#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 11:14:13 2019

@author: dobri
"""

import numpy as np
from astropy.stats import circmean

x = np.multiply(np.pi,[(0,1/4,2/4,3/4,4/4),(1,5/4,6/4,7/4,8/4),(5/4,5/4,5/4,5/4,5/4),(0/5,2/5,4/5,6/5,8/5)])
s = np.shape(x)

phikprime = np.array(x*0, dtype=complex);
phikprimebar = np.zeros((s[1],1), dtype=complex)
phikbar = np.zeros((s[0],1))
rhok = np.zeros((s[0],1))
    
for j in range(0,len(x)):
    for k in range(0,len(x[j,:])):
        phikprime[j,k]=np.complex(np.cos(x[j,k]),np.sin(x[j,k]))
    
    phikprimebar[j] = np.sum(phikprime[j,:])/s[1]
    phikbar[j] = np.angle(phikprimebar[j])
    rhok[j] = np.absolute(phikprimebar[j])
    print(phikbar[j],circmean(x[j,:]),rhok[j])