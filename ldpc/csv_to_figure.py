import numpy as np 
import csv
import matplotlib.pyplot as plt
from pylab import * 

R=np.zeros(7)
BER_amp = np.zeros(7)
BER_ldpc = np.zeros(7)
BER_ldpc_amp = np.zeros(7)
BER_sparc = np.zeros(7)

i=0
with open('RVsBER_amp_ldpc_3.csv') as myfile:
		reader = csv.DictReader(myfile)
		for row in reader:
			R[i] = row['R']
			BER_amp[i] = row['BER_amp']
			BER_ldpc[i] = row['BER_ldpc']
			BER_ldpc_amp[i] = row['BER_ldpc_amp']
			BER_sparc[i] = row['BER_sparc']
			i=i+1

print(R)
print(BER_amp)

P=2.4
sigma=1

capacity = 1/2*np.log2(1 + P/sigma**2)


fig, ax = plt.subplots()
ax.set_yscale('log', basey=10)
ax.plot(R, BER_amp, 'k--', label = 'SPARC with outer code: AMP only')
ax.plot(R, BER_ldpc, 'k:', label = 'SPARC with outer code: after LDPC')
ax.plot(R, BER_ldpc_amp, 'k-', label = 'SPARC with outer code: after final AMP')
ax.plot(R, BER_sparc, 'k-.', color='g', label = 'SPARC with no outer code')
plt.axvline(x=capacity, color='r', linestyle='-', label='Shannon limit')
plt.xlabel('Overall Rate')
plt.ylabel('BER')
plt.legend()
plt.savefig('RVsBER_amp_ldpc_3_fromcsv.png')