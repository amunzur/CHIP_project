#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 14 16:14:29 2021

@author: amunzur
"""
import pandas as pd
import numpy as np 

#################
# DECLARE VARIABLES
#################

PATH_both = "/groups/wyattgrp/users/amunzur/chip_project/data/GU-17-150-07Nov2017-cfDNA_WBC_filter_removed_somatic_anno_snv_EC.tsv"
PATH_tumor = "/groups/wyattgrp/users/amunzur/chip_project/data/GU-17-150-07Nov2017-cfDNA_somatic_anno_snv_EC.tsv"
PATH_wbc = "/groups/wyattgrp/users/amunzur/chip_project/data/wbc_variants.tsv"


both = pd.read_csv(PATH_both, sep = "\t")
tumor = pd.read_csv(PATH_tumor, sep = "\t")

# go through the both df and identify common lines
s1 = pd.merge(both, tumor, how='inner', on=['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene'])

# concat the dfs 
combined = pd.concat([both, tumor])
wbc = combined.drop_duplicates(subset = ['Chr', 'Start', 'End'], keep = False)

# save 
wbc.to_csv(PATH_wbc, sep = "\t")
