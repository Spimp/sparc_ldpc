# csv to figure code for soft_hard_plot()
import numpy as np 
import csv
import matplotlib.pyplot as plt
from pylab import * 
import re

# allows you pick out numbers from strings. 
numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
rx = re.compile(numeric_const_pattern, re.VERBOSE)

# Parameters you need to fill in about the simulation
datapoints = 10
soft_iter = 2
soft = True
hard = True
L=768
M=512
sec=569
logm=log(M)
r_sparc=1

EbN0_dB=np.zeros(datapoints)
BER_sparc = np.zeros(datapoints)
if soft:
	BER_ldpc_soft = np.zeros((datapoints, soft_iter))
	BER_amp_soft = np.zeros((datapoints, soft_iter+1))
if hard:
	BER_ldpc_hard = np.zeros(datapoints)
	BER_amp_hard = np.zeros((datapoints, 2))


i=0
j=0
with open('EbN0VsBER_soft_hard_100.csv') as myfile:
		reader_soft = csv.DictReader(myfile)
		for row in reader_soft:
			EbN0_dB[i] = row['EbN0_dB']
			BER_sparc[i] = row['BER_sparc']

			if soft: 
				BER_ldpc_soft[i,:] = rx.findall(row['BER_ldpc_soft'])
				BER_amp_soft[i,:] = rx.findall(row['BER_amp_soft'])
			if(i==datapoints-1):
				break
			i=i+1

		reader_hard = csv.DictReader(myfile)
		for row in reader_hard:
			if hard:	
				BER_ldpc_hard[j] = row['BER_ldpc_hard']
				BER_amp_hard[j,:] = rx.findall(row['BER_amp_hard'])
			j=j+1

# calculating the capacity so this can be plotted
nl=logm * sec
n = L*logm/r_sparc
R = (L*logm-nl*(1-5/6))/n
snrc = (2**(2*R) - 1)
EbN0c = 1/(2*R) * snrc
EbN0c_dB = 20*log10(EbN0c)



# plot the figure. Comment and uncomment lines as desired to plot the lines you want
fig, ax = plt.subplots()
ax.set_yscale('log', basey=10)
ax.plot(EbN0_dB, BER_sparc, 'g:', label = 'SPARC w/ same overall rate')
if soft:
    # comment out and add lines below as appropriate
    #ax.plot(EbN0_dB, BER_amp_soft[:,0], 'b:', label = 'SPARC w/ outer code: after AMP')
    #ax.plot(EbN0_dB, BER_ldpc_soft[:,0], 'k:', label = 'SPARC w/ outer code: after LDPC (soft)')
    #ax.plot(EbN0_dB, BER_amp_soft[:,1], 'b--', label = 'SPARC w/ outer code: after 2nd round of AMP (soft)')
    ax.plot(EbN0_dB, BER_ldpc_soft[:,1], 'k--', label = 'SPARC w/ outer code: after 2nd round of LDPC (soft)')
    #ax.plot(EbN0_dB, BER_amp_soft[:,2], 'b-', label = 'SPARC w/ outer code: after 3rd round of AMP (soft)')
if hard:
    #ax.plot(EbN0_dB, BER_amp_hard[:,0], 'c:', label = 'SPARC w/ outer code: after AMP (hard)')
    #ax.plot(EbN0_dB, BER_ldpc_hard, 'm:', label = 'SPARC w/ outer code: after LDPC (hard)')
    ax.plot(EbN0_dB, BER_amp_hard[:,1], 'b--', label = 'SPARC w/ outer code: after 2nd round of AMP (hard)')
plt.axvline(x=EbN0c_dB, color='r', linestyle='-', label='Shannon limit')
plt.xlabel('$E_b/N_0$ (dB)', fontsize=15) # at some point need to work out how to write this so it outputs properly
plt.ylabel('BER', fontsize=15)
plt.tight_layout()
plt.legend(loc=1, prop={'size': 9})
plt.savefig('EbN0VsBER_soft_hard_100_refined.png')
