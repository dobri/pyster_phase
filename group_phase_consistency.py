"""
Python implementation of the Frank & Richardson (2010) method for group phase
synchronization. They called it cluster phase, but it's better to call it 
group phase consistency, multivariate phase consistency, etc.. There's no 
clustering involved; this is confusing. You simply zero-center each 
variable by removing its mean phase relative to the grand mean phase. It 
doesn't matter what's the grand mean phase before zero-centering, after this 
operation it will be zero.

Attention! This method is biased by number of variables and data length. Don't
use it to compare b/w groups of different N or very different trial lengths.
See figures in the bias folder.

For significance testing of phase sync, as is always the case with these things,
don't just interpret the raw value of rho. You need to construct an 
appropriate null distribution using matching simulations or a surrogate procedure.

@author: dobri @LIVELab, McMaster University, 2019
@author: Victoria Cheesman
"""
import numpy as np
import csv
import sys, getopt
#from astropy.stats import rayleightest

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

if __name__ == '__main__':

    # Prepare the input and output file names and plotting flags from args.
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv,"hi:o:p:")
    except getopt.GetoptError:
        print('group_phase_consistency.py --help')
        sys.exit(2)
        
    if len(opts)==0:
        print('Usage: group_phase_consistency.py -i checkitworksdata4 -o <1 or 0> -p <1 or 0>')
        sys.exit(2)    


    LOGGING = False
    PLOTTING = False
    for opt, arg in opts:
        if opt == '-h':
            print('Usage: group_phase_consistency.py -i checkitworksdata4 -o <1 or 0> -p <1 or 0>')
            sys.exit()
        if opt in ("-i"):
            INPUT_FILE_NAME = arg
        if opt in ("-o"):
            LOGGING = bool(float(arg))
            OUTPUT_FILENAME = 'rhog_' + INPUT_FILE_NAME
        if opt in ("-p"):
            PLOTTING = bool(float(arg))


    # Import the data from a comma-delimited text file.
    data = []
    with open(INPUT_FILE_NAME) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            data.append(np.array(row,dtype='float'))
    data = np.array(data)
    s = np.shape(data)
    print('Imported data with %.0f' %s[0], 'rows (variables) by %.0f' %s[1], 'columns (observations).')
    
    # The core of the computation.
    qt, rt = circ_mean(data) # The group phase and vector length averaged per time step.
    phikt = data - np.transpose(qt) # Individual phases relative to the group, per time step.
    phikbar,rhok = circ_mean(np.transpose(phikt)) # The time-average of these individual to group relative phases.
    centered_phases = phikt - phikbar # Zero-center individuals' relative phases by subtracting their mean deviation from the group.
    qgroupt, rhogroupt = circ_mean(centered_phases) # The group phase and vector length averaged per time step from the centered rel phases.
    rhogroup = np.mean(rhogroupt) # Finally, the grand mean vector length.
    print('The group phase consistency is %.2f.' %rhogroup)
    
    # Save to a text file.
    if LOGGING:
        f = open(OUTPUT_FILENAME,'w')
        f.write("%.3f,\n" %rhogroup)
        f.close()
    
    # Plotting.
    if PLOTTING:
        import matplotlib.pyplot as plt
        plt.plot(rhogroupt,'-k')
        plt.ylabel(r'$ \rho _{group,i}$')
        plt.show()
