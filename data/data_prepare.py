#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
    File Name: data_prepare.py
    Description:
    
Created by YongBai on 2020/12/2 5:04 PM.
"""

import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import logging


class STSimpleData:

    def __init__(self, adata):

        self.adata = adata

    @classmethod
    def get_simple_st(cls, libary_id):

        # load data from scanpy
        adata = sc.datasets.visium_sge(libary_id)
        adata.var_names_make_unique()
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        adata.var['mt'] = [gene.startswith('MT-') for gene in adata.var_names]
        adata.obs['mt_frac'] = adata[:, adata.var['mt']].X.sum(1).A.squeeze() / adata.obs['total_counts']

        # QC and data preprocess
        sc.pp.filter_cells(adata, min_counts=5000)
        logging.info(f'Number of cells after min count filter: {adata.n_obs}')
        sc.pp.filter_cells(adata, max_counts=35000)
        logging.info(f'Number of cells after max count filter: {adata.n_obs}')
        adata = adata[adata.obs['mt_frac'] < 0.2]
        logging.info(f'Number of cells after MT filter: {adata.n_obs}')
        sc.pp.filter_cells(adata, min_genes=3000)
        logging.info(f'Number of cells after gene filter: {adata.n_obs}')
        sc.pp.filter_genes(adata, min_cells=10)
        logging.info(f'Number of genes after cell filter: {adata.n_vars}')

        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=2000, inplace=True)

        sc.pp.pca(adata)

        return cls(adata)

