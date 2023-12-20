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

When computing nodal fraction for hydrogen and ammonia, missing values are filled using fractions from electricity.

Data are in TWh in this rule.
"""

import logging

import pandas as pd
import numpy as np

from _helpers import configure_logging

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industrial_energy_demand_patex",
                                   simpl="",
                                   clusters="47",
                                   planning_horizons=2030)

    configure_logging(snakemake)

    industrial_demand = pd.read_csv(snakemake.input["industrial_demand"], index_col=[0])
    industry_patex = pd.read_csv(snakemake.input["industry_patex"], index_col=[0]).T
    industry_patex.index = pd.MultiIndex.from_tuples(
        industry_patex.index.str.split("_", n=1).to_list(), names=["country", "sector"])
    industry_patex = industry_patex.reset_index().pivot_table(
        values=snakemake.wildcards.planning_horizons, index="country", columns="sector")

    industry_patex["share_methanol"] = (
            industry_patex["methanol"]
            / (industry_patex["methanol"] + industry_patex["fischer_tropsch"])
    ).fillna(0)

    industry_patex["electricity"] = industry_patex["elec_ind"]
    industry_patex["hydrogen"] = (
            industry_patex["h2_sectors"]  # demand from sectors
            - industry_patex["h2_ind_nh3"]  # used to produce ammonia
            + industry_patex["h2_efuels"] * (1 - industry_patex["share_methanol"])  # used for efuels but shipping
    )
    industry_patex["ammonia"] = (
            industry_patex["nh3_ind"] * snakemake.config["industry"]["MWh_NH3_per_tNH3"] * 1e-3  # industrial demand
            + (industry_patex["nh3_marine"]  # shipping demand
               * snakemake.config["sector"]["shipping_methanol_share"][snakemake.wildcards.planning_horizons])
    )

    carriers = ["electricity", "hydrogen", "ammonia"]
    carriers_excluded = [c for c in industrial_demand.columns if c not in carriers]
    industry_patex = industry_patex[carriers]
    industry_patex[carriers_excluded] = 0

    industrial_demand["country"] = industrial_demand.index.str[:2]
    industrial_demand.set_index("country", append=True, inplace=True)
    nodal_fraction = (
            industrial_demand
            / industrial_demand[[]].join(industrial_demand.groupby(by="country").sum(), on="country")
    )
    nodal_fraction = nodal_fraction.replace([np.inf, -np.inf], np.nan).fillna({
        "hydrogen": nodal_fraction["electricity"],
        "ammonia": nodal_fraction["electricity"],
    })
    nodal_fraction[carriers_excluded] = 0

    industrial_demand_patex = (
            industrial_demand[[]].join(industry_patex, on="country")
            * nodal_fraction
    ).droplevel(1)

    [logging.info(
        f"Industrial {c} load is {v:.2f} TWh") for c, v in industrial_demand_patex[carriers].sum().to_dict().items()]
    industrial_demand_patex.to_csv(snakemake.output["industrial_demand_patex"])
