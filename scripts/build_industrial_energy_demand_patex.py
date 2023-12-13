# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2023 Climact for The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

"""
This rule override the industrial energy demand defined by default sector module using data from the 2050
Pathways Explorer. Data used are to override are the electrical, ammonia and hydrogen industrial energy demand.

In the 2050 Pathways Explorer, ammonia for shipping is not modeled. The energy demand for ammonia is replaced by a
demand in methanol. To overcome this limitation in PyPSA (where ammonia is well modeled), we shift the demand for H2
used for methanol in ammonia demand. PyPSA will then choose the best way to produce this ammonia.
"""

import logging

import pandas as pd

from _helpers import configure_logging

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industrial_energy_demand_patex",
                                   simpl="",
                                   clusters="47",
                                   planning_horizons=2030)

    configure_logging(snakemake)

    years = snakemake.config["scenario"]["planning_horizons"]
    index_h2 = [str(y) + "_hydrogen" for y in years]
    index_nh3 = [str(y) + "_ammonia" for y in years]
    index_elec = [str(y) + "_electricity" for y in years]
    industry = df_results[["region", "sector"] + years]
    countries = get_eu27_countries()

    result = pd.DataFrame([[] * 3], index=index_h2 + index_nh3 + index_elec)
    for co in countries:
        df_co = industry[industry.region.isin([co])].drop(columns=["region"]).set_index("sector")
        share_e_methanol = (df_co.loc["methanol"] / (df_co.loc["methanol"] + df_co.loc["fischer_tropsch"])).fillna(0)
        h2_demand = df_co.loc["h2_sectors"] - df_co.loc["h2_ind_nh3"] + df_co.loc["h2_efuels"] * (1 - share_e_methanol)
        nh3_demand = df_co.loc["nh3_ind"] * 5.166e-3 + df_co.loc["nh3_marine"] * \
                     pd.Series(snakemake.config["sector"]["shipping_methanol_share"]).loc[years]
        elec_demand = df_co.loc["elec_ind"]
        result.loc[index_h2, co] = h2_demand.values
        result.loc[index_nh3, co] = nh3_demand.values
        result.loc[index_elec, co] = elec_demand.values

    df_results = result

    industry = pd.read_csv(snakemake.input["industrial_demand"], index_col=[0])
    patex = pd.read_csv(snakemake.input["industry_patex"], index_col=[0]).T
    new_demand = industry * 0
    for carrier in ["hydrogen", "ammonia", "electricity"]:
        for country in industry.index.str[:2].unique():
            ca_co = industry.loc[new_demand.filter(like=country, axis=0).index, carrier]
            if abs(ca_co.sum()) > 1e-6:
                weights = (ca_co / ca_co.sum()).fillna(0)
            else:
                ca_co = industry.loc[new_demand.filter(like=country, axis=0).index, 'electricity']
                weights = (ca_co / ca_co.sum()).fillna(0)

            patex_co = patex.loc[patex.filter(like=country, axis=0).index, carrier].values
            new_demand.loc[new_demand.filter(like=country, axis=0).index, carrier] = weights * patex_co
        logging.info(f"Carrier {carrier} amounts to {sum(new_demand[carrier])} [TWh]")
    new_demand.to_csv(snakemake.output["industry_demand"])
