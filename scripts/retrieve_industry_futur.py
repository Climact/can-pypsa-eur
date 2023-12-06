# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
This rule downloads yearly electric load data for each european country from the `2050 Pathways Explorer <https://pathwaysexplorer.climact.com>`_.
This rule downloads appropriate years using configuration file. Yearly electric load data are defined by per-country
scenarios.

**Releveant Settings**

.. code:: yaml

    planning_horizons

    snapshots:
        start:

**Inputs**

- ``scenario_builder_tool_input.xlsx`` Projection scenario definitions for each region of Europe.
- ``resources/load.csv`` Hourly per-country load profiles.

**Outputs**

- ``data/patex`` Yearly electric load per-country.
"""

import json
import logging
import os
import re
from pathlib import Path

import pandas as pd
import requests

from _helpers import configure_logging
from retrieve_load_futur import (__safely_to_int,
                                 load_scenario_builder,
                                 get_eu27_countries,
                                 fill_scenario_list,
                                 get_results)



# ToDo Move this to yaml configuration
URL_TO_USE = "prod url"
# URL_TO_USE = "test url"

METRIC_MAP = pd.DataFrame([
    ["clm_CO2-need-by-way-of-prod_CCU_elc_e-methanol[MtCO2e]", "methanol"],
    ["clm_CO2-need-by-way-of-prod_CCU_elc_fischer-tropsch[MtCO2e]","fischer_tropsch"],
    ["elc_energy-demand-by-direct-use-and-energy-carrier_for-sector_hydrogen[TWh]","h2_sectors"],
    ["ind_energy-demand-by-carrier-feedstock-group_hydrogen_feedstock_chemicals-group[TWh]","h2_ind_nh3"],
    ["elc_energy-demand-by-direct-use-and-energy-carrier_for-power-prod_hydrogen[TWh]","h2_efuels"],
    ["elc_energy-demand-by-direct-use-and-energy-carrier_for-power-prod_hydrogen[TWh]","h2_methanol"],
    ["ind_material-production-by-material_chemical-ammonia[kt]","nh3_ind"],
    ["tra_energy-demand-bunkers-by-type_bunkers-marine[TWh]","nh3_marine"],
    ["ind_energy-demand-by-carrier_electricity[TWh]", "elec_ind"],
], columns=["metric_id", "sector"])

TEMPLATE_REQUEST = """
{
  "levers":
  [LEVERS],
  "outputs": {"0": [VARIABLES]},
  "dimension": {"0": null},
  "aggregate": true
}
"""
API = {
    "prod url": {
        "db": "https://pathwaysexplorer.climact.com/api/v1.0/model_results",
        "model": "https://pathwaysexplorer.climact.com/model/v1.0/model_results"
    },
    "test url": {
        "db": "https://test.pathwaysexplorer.climact.com/api/v1.0/model_results",
        "model": "https://test.pathwaysexplorer.climact.com/model/v1.0/model_results"
    }
}

RX_SCENARIO = re.compile("^\[([A-z]+)\]\s(.*)$")

def format_results(results):
    """
    Format received JSON into dataframe
    Change unit from TWh to MWh (* 1e6)
    Add non-MS data using proportions

    :param results: data in JSON format
    :return df_results: data in dataframe
    """
    df_results = []
    for (master_region, scenario), values in results.items():
        values = values['outputs']['0']
        timeAxis = values[0]['timeAxis']
        data = []
        names = []
        for i in values:
            data.append(i['data'])
            if RX_SCENARIO.match(scenario):
                region = RX_SCENARIO.match(scenario).group(1)
            else:
                region = master_region
            names.append((scenario, region, i['id'], i['title']))
        df_results.append(
            pd.DataFrame(data, index=pd.MultiIndex.from_tuples(names, names=['scenario', 'region', 'metric_id',
                                                                             'metric_name']),
                         columns=timeAxis))

    df_results = pd.concat(df_results).reset_index(drop=False)

    df_results = (
        df_results.set_index("metric_id")
        .join(METRIC_MAP.set_index("metric_id"))
        .set_index("sector").reset_index()
    )
    # Hypothesis : One unique scenario for each country
    df_results = df_results.groupby(by=["region", "sector"]).sum().reset_index()
    
    years = snakemake.config["scenario"]["planning_horizons"]
    index_h2 = [str(y) +"_hydrogen" for y in years]
    index_nh3 = [str(y) +"_ammonia" for y in years]
    index_elec = [str(y) +"_electricity" for y in years]
    industry = df_results[["region","sector"]+years]
    countries = get_eu27_countries
    
    
    result = pd.DataFrame([[]*3],index = index_h2+index_nh3+index_elec)
    for co in countries:
        df_co = industry[industry.region.isin([co])].drop(columns=["region"]).set_index("sector")
        share_e_methanol = (df_co.loc["methanol"] / (df_co.loc["methanol"] + df_co.loc["fischer_tropsch"])).fillna(0)
        h2_demand =  df_co.loc["h2_sectors"] - df_co.loc["h2_ind_nh3"] + df_co.loc["h2_efuels"] - share_e_methanol * df_co.loc["h2_methanol"]
        nh3_demand = (df_co.loc["nh3_ind"]*5.166e-3 + df_co.loc["nh3_marine"])* [0,0.68,1]
        elec_demand = df_co.loc["elec_ind"]
        result.loc[index_h2,co] = h2_demand.values
        result.loc[index_nh3, co] = nh3_demand.values
        result.loc[index_elec, co] = elec_demand.values
    
    df_results = result
    df_results.columns = df_results.columns.str.replace("EL", "GR")

    # Manual scaling for non-MS if missing
    dict_non_MS = {"AL": "GR",
                   "MK": "GR",
                   "NO": "SK",
                   "CH": "FR",
                   "GB": "FR",
                   "BA": "HR",
                   "KV": "HU",
                   "RS": "HU",
                   "ME": "GR"}
    dict_non_MS = {k: v for k, v in dict_non_MS.items() if not any(df_results.columns.str.startswith(k))}
    historical_load_h = pd.read_csv(snakemake.input.load_hourly, index_col="utc_timestamp", parse_dates=True)
    missing_load = []
    
    for non_MS, MS in dict_non_MS.items():
        df = df_results.loc[:,df_results.columns.str.startswith(MS)]
        df *= historical_load_h[non_MS].sum() / historical_load_h[MS].sum()
        df.columns = df.columns.str.replace(MS, non_MS)
        missing_load.append(df)
    df_results = pd.concat([df_results] + missing_load,axis=1)

    return df_results


def write_files(df_results):
    """
    Export of data
    """
    df_results.index = df_results.index.set_names("year")

    path = Path(snakemake.output[0]).parent
    for y in snakemake.config["scenario"]["planning_horizons"]:
        try:
            df_to_write = df_results.filter(like=str(y),axis=0)
            df_to_write.index = df_to_write.index.str.replace(f"{y}_","")
            df_to_write.to_csv(Path(path, f"patex_ind_{y}.csv"))
        except OSError:
            os.makedirs(path)
            df_results.loc[[y]].to_csv(Path(path, f"patex_ind_{y}.csv"))


if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_industry_futur")

    configure_logging(snakemake)

    # Load configuration
    logging.info("Loading configuration")
    df_scenarios = load_scenario_builder()
    scenarios_dict = fill_scenario_list(df_scenarios)

    # Getting data from API
    logging.info("Getting data from API")
    results = get_results(scenarios_dict, METRIC_MAP["metric_id"])

    # Formatting data
    logging.info("Formatting data")
    df_results = format_results(results)

    # Writing data
    logging.info("Writing data")
    write_files(df_results)