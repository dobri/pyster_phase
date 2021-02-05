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
#from astropy.stats import circmean
import matplotlib.pyplot as plt
import datetime


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


def cluster_phase(data):
    qt, rt = circ_mean(data) # The group phase and vector length averaged per time step.
    phikt = data - np.transpose(qt) # Individual phases relative to the group, per time step.
    phikbar,rhok = circ_mean(np.transpose(phikt)) # The time-average of these individual to group relative phases.
    centered_phases = phikt - phikbar # Zero-center individuals' relative phases by subtracting their mean deviation from the group.
    centered_phases = np.mod(centered_phases,np.pi*2)
    # There's no point of this. By centering each variable we created a very strong central bias.
    #from astropy.stats import rayleightest
    #s = np.shape(data)
    #rp = rayleightest(np.reshape(centered_phases,(s[0]*s[1],1)))
    #print('The Rayleigh test is p=%.4f' %rp)
    qgroupt, rhogroupt = circ_mean(centered_phases) # The group phase and vector length averaged per time step from the centered rel phases.
    rhogroup = np.mean(rhogroupt) # Finally, the grand mean vector length.
    return rhogroup

def sine_with_rand_phases_periods(N,n,fr=1,var_fr=0):
    dt = .01
    t = np.multiply(range(1,n+1),dt)
    data = np.zeros([len(t),N])
    for c in range(0,N):
        phase = np.random.rand()*2*np.pi
        fr = fr + var_fr*(np.random.rand()*2-1)
        tt = t*fr*2*np.pi+phase + np.random.randn(1,len(t))*.5
        data[:,c] = np.angle(1j*np.sin(tt)+1*np.cos(tt))
    return data


if __name__ == '__main__':
    # 2^p time points, nump participants.
    #powers = np.linspace(3,12,6,dtype='int16')
    powers = np.linspace(8,14,7,dtype='int16') #7
    num_participants = np.power(2,range(1,5)) # 8
    cent_freq = (1,2,4,10)
    sine_or_random = 1
    randfreq = range(0,2)
    for ran in randfreq:
        for nump in num_participants:
            fig = plt.figure()
            pp_loop = 0
            for fr in cent_freq:
                rhogroup_vec_sum1 = np.zeros_like(powers,dtype='float')
                rhogroup_vec_sum2 = np.zeros_like(powers,dtype='float')
                for rep in range(1,101): # repetitions for smoothness and stats
                    c=0
                    for p in powers:
                        if (sine_or_random == 1):
                            data = sine_with_rand_phases_periods(nump,2**p,fr,ran)
                        else:
                            data = 2*np.pi*np.random.rand(nump,2**p)
                        #s = np.shape(data)
                        #print('Imported data with %.0f' %s[0], 'rows (variables) by %.0f' %s[1], 'columns (observations).')
                        
                        # The core of the computation.
                        rhogroup = cluster_phase(data)
                        rhogroup_vec_sum1[c] = rhogroup_vec_sum1[c] + rhogroup
                        rhogroup_vec_sum2[c] = rhogroup_vec_sum2[c] + rhogroup**2
                        c = c+1
                        print('N=%2.0f' %nump,', n=%8.0f' %(2**p), ', rep=%4.0f' %rep, ', cpc=%4.2f' %rhogroup)
        
                sem = ((rhogroup_vec_sum2 - rhogroup_vec_sum1**2/rep)/(rep-1))**.5/(rep**.5)
                m = rhogroup_vec_sum1/rep # The average rho across repetitions.
                
                pp_loop = pp_loop + 1
                plt.errorbar(powers+(np.random.rand()-.5)*.2, m, yerr=sem, label=fr, color=plt.cm.copper(pp_loop / len(cent_freq)))
                #plt.plot(powers, m, '-ok', color=plt.cm.copper(pp_loop / len(num_participants)))
        
            plt.xlabel('$2^n$ number of time points in the series.')
            plt.ylabel(r'$\rho_{group}$')
            plt.ylim((0,1))
            plt.legend(cent_freq)
            filenamestring = 'rhopic_rand'+str(ran)+'_'+str(nump)+'pps'
            plt.title(filenamestring)
            plt.draw()
            filenamestring = filenamestring+'_'+datetime.datetime.now().strftime("%Y%B%d-%I%M%S%p")+'.png'
            plt.savefig(filenamestring, dpi=300)
            #plt.show()

    """
    plt.plot(np.log2((4,6,8,10,16,32,64,128,256)),(.45,.36,.32,.29,.22,.16,.11,.08,.06),'-ok')
    plt.xlabel('$2^n$ number of participants in the group')
    plt.ylabel(r'$\rho_{group}$')
    plt.title('Cluster phase from random (uniform) series')
    """
    
    """
    What's going on here? There's a serious mistake somehwere. The polar 
    histograms seem strongly warped to the center. But, importantly, I get the
    same rho_i (over time) as they have in their paper (~.4) in the uncoupled
    condition. Does this mean that they've been reporting this flawed analysis
    all along? Scary!
    
    It seems that taking relative phase in each time step while having a low
    number of variables (participants) creates the possibility that rel phases
    are always clustered. As you increase the variables n, the histogram doesn't
    look so biased any more, but it's still there.
    
    Re-centering inflates consistency. Even if there's a small central tendency 
    in each variable these will be aligned.
    """