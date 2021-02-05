"""
Python implementation of the Frank & Richardson (2010) method for group phase
synchronization. They called it cluster phase, but it's better to just call it 
group phase consistency, or multivariate phase consistency. There's no 
clustering involved, this is a confusing term. We simply zero-center each 
variable by removing its mean phase relative to the grand mean phase. It 
doesn't matter what's the grand mean phase before zero-centering, after this 
operation it will be zero.

@author: dobri @LIVELab, McMaster University, 2019
"""
import numpy as np
#import csv
#import sys, getopt
from astropy.stats import rayleightest
from astropy.stats import circmean
import matplotlib.pyplot as plt


def circ_mean(x):
    s = np.shape(x)
    aprime = np.array(x*0, dtype=complex);
    a = np.zeros((s[1],1), dtype=float)
    r = np.zeros((s[1],1), dtype=float)
    for i in range(0,s[1]):
        for j in range(0,s[0]):
            aprime[j,i]=np.complex(np.cos(x[j,i]),np.sin(x[j,i]))
        aprimebar = np.sum(aprime[:,i])/s[0]
        a[i] = np.angle(aprimebar)
        r[i] = np.absolute(aprimebar)
    return a, r


def polar_histogram(ang,s):
    bin_n = 100
    counts = np.zeros((s[0],bin_n), dtype=float)
    bins = np.zeros((s[0],bin_n), dtype=float)
    i=0
    for d in ang:
        c,b=np.histogram(d,bins = bin_n)
        counts[i,] = c
        bins[i,] = b[0:len(b)-1]
        i=i+1

    bottom = np.min(counts)
    width = (2*np.pi) / bin_n

    fig = plt.figure()
    ax = fig.gca(polar=True)
    for i in range(0,s[0]):
        bars = ax.bar(bins[i,], counts[i,], width=width, bottom=bottom)
        bars[i].set_facecolor(plt.cm.jet(i / s[0]))
        bars[i].set_alpha(0.2)
        print(i)

    plt.show()
    return


if __name__ == '__main__':
    """
    This simple centering solves the warped histogram issue. Is this a smarter
    way of measuring group synchrony?
    """

    data = 2*np.pi*np.random.rand(6,1000)
    s = np.shape(data)
    print('Imported data with %.0f' %s[0], 'rows (variables) by %.0f' %s[1], 'columns (observations).')

    data_cent = data*0
    i = 0
    for d in data:
        d = np.mod(d - circmean(d),np.pi*2)
        data_cent[i,] = d
        i = i+1
    
    qgroupt, rhogroupt = circ_mean(data_cent)
    plt.plot(rhogroupt,'-k')
    plt.ylabel(r'$ \rho _{group,i}$')
    plt.show()

    polar_histogram(data,s)
    rp0 = rayleightest(np.reshape(data,(s[0]*s[1],1)))
    print('The Rayleigh test is p=%.4f.' %rp0)

    polar_histogram(data_cent,s)
    rp1 = rayleightest(np.reshape(data_cent,(s[0]*s[1],1)))
    print('The Rayleigh test is p=%.4f' %rp1)

    for d1,d2 in zip(data,data_cent):
        rp2 = rayleightest(d1)
        rp3 = rayleightest(d2)
        print('The Rayleigh tests are p=%.4f and %.4f.' %(rp2,rp3))
        
