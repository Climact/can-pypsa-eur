# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
This rule downloads yearly industrial load data for each european country from the
`2050 Pathways Explorer <https://pathwaysexplorer.climact.com>`_. This rule consider as industrial load: electricity,
ammonia and hydrogen. This rule downloads appropriate years using configuration file. Yearly load data are defined
by per-country scenarios.

This rule mainly reuse definitions from ``retrieve_load_futur``.

**Releveant Settings**

.. code:: yaml

    planning_horizons

    snapshots:
        start:

    snapshots:
        start:

**Inputs**

- ``scenario_builder_tool_input.xlsx`` Projection scenario definitions for each region of Europe.
- ``resources/load.csv`` Hourly per-country load profiles.

**Outputs**

- ``data/patex_ind`` Yearly industrial load per-country.
"""

import logging

import pandas as pd

from _helpers import configure_logging
from retrieve_load_futur import fill_scenario_list
from retrieve_load_futur import format_results
from retrieve_load_futur import get_results
from retrieve_load_futur import load_scenario_builder
from retrieve_load_futur import write_files

METRIC_MAP = pd.DataFrame([
    ["clm_CO2-need-by-way-of-prod_CCU_elc_e-methanol[MtCO2e]", "methanol"],
    ["clm_CO2-need-by-way-of-prod_CCU_elc_fischer-tropsch[MtCO2e]", "fischer_tropsch"],
    ["elc_energy-demand-by-direct-use-and-energy-carrier_for-sector_hydrogen[TWh]", "h2_sectors"],
    ["ind_energy-demand-by-carrier-feedstock-group_hydrogen_feedstock_chemicals-group[TWh]", "h2_ind_nh3"],
    ["elc_energy-demand-by-direct-use-and-energy-carrier_for-power-prod_hydrogen[TWh]", "h2_efuels"],
    ["ind_material-production-by-material_chemical-ammonia[kt]", "nh3_ind"],
    ["tra_energy-demand-bunkers-by-type_bunkers-marine[TWh]", "nh3_marine"],
    ["ind_energy-demand-by-carrier_electricity[TWh]", "elec_ind"],
], columns=["metric_id", "sector"])

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("retrieve_industry_futur",
                                   configfiles="config.CANEurope.runner.yaml")

    configure_logging(snakemake)

    # Load configuration
    logging.info("Loading configuration")
    df_scenarios = load_scenario_builder(snakemake.input.scenario_builder)
    scenarios_dict = fill_scenario_list(df_scenarios)

    # Getting data from API
    logging.info("Getting data from API")
    results = get_results(scenarios_dict, METRIC_MAP["metric_id"])

    # Formatting data
    logging.info("Formatting data")
    horizons = [pd.Timestamp(snakemake.config["snapshots"]["start"]).year] + \
               snakemake.config["scenario"]["planning_horizons"]
    df_results = format_results(results, horizons, METRIC_MAP, snakemake.input.load_hourly) * 1e-6  # TWh

    # Writing data
    logging.info("Writing data")
    write_files(df_results, snakemake.output, "patex_ind_")
