#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 23 03:16:08 2019

@author: dobri
"""

import numpy as np
from astropy.stats import circmean
from astropy import units as u
import csv

data = []
with open('data.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        data.append(np.array(row,dtype='float'))
        print(data[line_count])
        print(circmean(data[line_count]/360*2*np.pi)*u.rad)
        print(circmean(data[line_count]/360*2*np.pi)/np.pi/2*360*u.deg)
        line_count += 1
    
    print(f'Processed {line_count} lines.')