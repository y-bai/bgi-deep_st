#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
    File Name: train.py
    Description:
    
Created by YongBai on 2020/10/14 2:21 PM.
"""

import os
import argparse
import torch
import torch.nn as nn
import logging
import scanpy as sc
from logger import get_logger
from prepare import *
from data import *


from torchsummary import summary

# fix random seeds for reproducibility
SEED = 123
torch.manual_seed(SEED)  # CPU seed
torch.cuda.manual_seed(SEED)  # GPU seed
# accelerate pytorch
# using deterministic cudnn convolutional opt(like seed = 0 so that all results of cnn can be reproducible)
torch.backends.cudnn.deterministic = True
# if input and model on the computation graph not changing,
# then benchmark = True would help speed up pytorch
torch.backends.cudnn.benchmark = True


def main(in_args):

    # setup paramters
    config = Config.init_params(in_args.config).config

    # get logger
    logger = get_logger(__name__)

    logger.debug(config)
    logger.info('Pytorch version: {}'.format(torch.__version__))
    logger.info('cuda version:{}'.format(torch.version.cuda))

    # adata = STSimpleData.get_simple_st('V1_Human_Lymph_Node').adata

    # print(adata)
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()

if __name__ == '__main__':
    args = argparse.ArgumentParser(description='SpatialTranscript')
    args.add_argument('-c', '--config',
                      default=None, type=str,
                      help='config json file path')
    args = args.parse_args()
    main(args)