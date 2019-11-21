#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 10:04:10 2019

@author: yk
"""

import os
import logging
import argparse

def mkdir(dir_name):
    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)
        logging.info(' mkdir '+ dir_name ) 
        
if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s',level=logging.INFO)
    parse = argparse.ArgumentParser(description='make a new dir')
    parse.add_argument('dirname')
    args = parse.parse_args()
    
    mkdir(args.dirname)
    