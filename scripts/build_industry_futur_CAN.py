# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 17:39:53 2023

@author: VincentLaguna
"""

import logging
import os
import re
from pathlib import Path

import pandas as pd
from _helpers import configure_logging

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_industry_futur_CAN",
                                   simpl = "",
                                   clusters="47",
                                   planning_horizons=2030)

    configure_logging(snakemake)
    
    
    industry = pd.read_csv(snakemake.input["industry_demand"],index_col=[0])
    patex = pd.read_csv(snakemake.input["patex_industry"],index_col=[0]).T
    new_demand = industry * 0
    for carrier in ["hydrogen","ammonia","electricity"]:
        for country in industry.index.str[:2].unique():
            ca_co = industry.loc[new_demand.filter(like=country,axis=0).index,carrier]
            weights = 0
            if abs(ca_co.sum())>1e-6:
                weights = (ca_co/ca_co.sum()).fillna(0)
            
            patex_co = patex.loc[patex.filter(like=country,axis=0).index,carrier].values
            new_demand.loc[new_demand.filter(like=country,axis=0).index,carrier] = weights*patex_co
            
    new_demand.to_csv(snakemake.output["new_demand"])
