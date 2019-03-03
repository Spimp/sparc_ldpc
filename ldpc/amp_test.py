import numpy as np
import math as math
from pylab import * 
import py.ldpc as ldpc
import matplotlib.pyplot as plt
import time
import csv
from bitarray import bitarray
from sparc_ldpc import * 

# Testing the amp function to ensure that it does give improvements in the error performance when given a non-zero beta_0

# the amp algorithm
def amp_test(y, σ_n, Pl, L, M, T, Ab, Az, β=np.array([None])):
    P = np.sum(Pl)
    n = y.size
    if β.all()==None:
        β = np.zeros((L*M, 1))
        z = y
        last_τ = 0
    else:
        β = β.reshape(L*M,1)
        z = y - Ab(β)
        # just to ensure the if statement isn't true on the first loop
        last_τ = 0
    # delete this line! Just for testing purposes
    #T = 10
    
    for t in range(T):
        τ = np.sqrt(np.sum(z**2)/n)
        if τ == last_τ:
        # can change params for isclose if need be. If stopping too early.
        #if np.isclose(τ, last_τ):
            #print("Stopping after t runs:", t)
            return β, t
        last_τ = τ

        
        # added astype to increase precision to avoid divide by zero in LLR
        s = (β + Az(z))#.astype(np.longdouble)
        rt_n_Pl = np.sqrt(n*Pl).repeat(M).reshape(-1, 1)
        u = s * rt_n_Pl / τ**2
        max_u = u.max()
        exps = np.exp(u - max_u)
        sums = exps.reshape(L, M).sum(axis=1).repeat(M).reshape(-1, 1)
        β = (rt_n_Pl * exps / sums).reshape(-1, 1)

        z = y - Ab(β) + (z/τ**2) * (P - np.sum(β**2) / n)
    
    return β, t

# comparing initialising the amp decoding with all zeros and initialising with the correct beta. 
def amp_init_test(L, M, snr_dB, P, r_sparc):
	logm = np.log2(M)
	total_bits = int(L*logm)
	snr = 10**(snr_dB/20)
	sigma = np.sqrt(P/snr)
	n = int(L*np.log2(M) / r_sparc)
	# uniform power allocation across sections
	Pl = P/L * np.ones(L)
	# generate beta

	# Generate the bits required
	sparc_bits = np.random.randint(0, 2, total_bits).tolist()

	assert len(sparc_bits)==total_bits

	# convert the bits to indices for encoding
	sparc_indices = bits2indices(sparc_bits, M)

	assert len(sparc_indices)==L
	   
	# Generate the SPARC transform functions A.beta and A'.z
	Ab, Az, ordering = sparc_transforms(L, M, n)

	# Generate our transmitted signal X
	beta = np.zeros((L*M, 1))
	for l in range(L):
		beta[l*M + sparc_indices[l]] = np.sqrt(n * Pl[l])
	x = Ab(beta)

	# transmit this x over the channel to produce y
	# Generate random channel noise and then received signal y
	z = np.random.randn(n, 1) * sigma
	y = (x + z).reshape(-1, 1)
	# distort beta in some small way and set this to beta_0
	beta_0 = beta 
	# run amp decoding of y with this beta_0 and also with an all zero beta_0 
	# output of amp with the initialisation beta_0
	beta_init, t_init = amp_test(y, 0, Pl, L, M, 64, Ab, Az, beta_0)
	beta_no_init, t_no_init = amp_test(y, 0, Pl, L, M, 64, Ab, Az)
	# compare the error performance and the number of iterations before termination. 

	# Convert decoded beta back to a message
	rx_message_init = []
	rx_message_no_init = []
	for l in range(L):
		idx = np.argmax(beta_init[l*M:(l+1)*M])
		rx_message_init.append(idx)
		idx2 = np.argmax(beta_no_init[l*M:(l+1)*M])
		rx_message_no_init.append(idx2)

	# Compute BER 
	ber_init = []
	ber_no_init = []
	ber_init.append(sum(bin(a^b).count('1') for (a, b) in zip(sparc_indices, rx_message_init))/total_bits)
	ber_no_init.append(sum(bin(a^b).count('1') for (a, b) in zip(sparc_indices, rx_message_no_init))/total_bits)
	print("For initialised amp, BER= ", ber_init, " and iterations= ", t_init)
	print("For amp with all zero beta_0, BER= ", ber_no_init, " and iterations= ", t_no_init)
	return ber_init, ber_no_init

def bp2sp_effect(L, M, r_sparc, P, sigma):
	# Look at what bp2sp is doing to the distribution of probabilities to see if I could correct for this. 
	# uniform power allocation
	Pl = P/L * np.ones(L)
	# start with some value of beta
	logm=np.log2(M)
	total_bits=int(L*logm)
	n = int(L*np.log2(M) / r_sparc)
	# Generate the bits required
	sparc_bits = np.random.randint(0, 2, total_bits).tolist()
	# convert the bits to indices for encoding
	sparc_indices = bits2indices(sparc_bits, M)
	# Generate the SPARC transform functions A.beta and A'.z
	Ab, Az, ordering = sparc_transforms(L, M, n)

	# Generate our transmitted signal X
	beta = np.zeros((L*M, 1))
	for l in range(L):
		beta[l*M + sparc_indices[l]] = np.sqrt(n * Pl[l])
	x = Ab(beta)

	# transmit this x over the channel to produce y
	# Generate random channel noise and then received signal y
	z = np.random.randn(n, 1) * sigma
	y = (x + z).reshape(-1, 1)

	beta_T, _= amp_test(y, 0, Pl, L, M, 64, Ab, Az)
	beta_T = beta_T/np.sqrt(n*P/L)
	# convert to bitwise posterior
	bp = sp2bp(beta_T, L, M)
	# convert back to sectionwise posterior and replot
	sp = bp2sp(bp, L, M)


	fig1, ax1 = plt.subplots()
	ax1.plot(beta_T, 'g--')
	
	fig, ax = plt.subplots()
	ax.plot(beta_T, 'g--', label='before')
	ax.plot(sp, 'b:', label='after')
	plt.ylabel("Sectionwise probability")
	plt.xlabel("Position")
	plt.legend()
	plt.show()





if __name__ == "__main__":
	L=512
	M=512
	# determines the number of sections which are set as zero in beta_0. The remaining sections are set as correctly decoded. 
	L_zero = 154
	P=4
	snr_dB=10
	snr = 10**(snr_dB/20)
	sigma = np.sqrt(P/snr)
	r_sparc = 1
	T=64
	repeats=100

	# uniform power allocation
	Pl = P/L * np.ones(L)
	# start with some value of beta
	logm=np.log2(M)
	total_bits=int(L*logm)
	n = int(L*np.log2(M) / r_sparc)
	ber_hard = 0
	ber_soft = 0
	ber_no_init = 0
	for i in range(repeats):
		# Generate the bits required
		sparc_bits = np.random.randint(0, 2, total_bits).tolist()
		# convert the bits to indices for encoding
		sparc_indices = bits2indices(sparc_bits, M)
		# Generate the SPARC transform functions A.beta and A'.z
		Ab, Az, ordering = sparc_transforms(L, M, n)

		# Generate our transmitted signal X
		beta = np.zeros((L*M, 1))
		for l in range(L):
			beta[l*M + sparc_indices[l]] = np.sqrt(n * Pl[l])
		x = Ab(beta)

		# transmit this x over the channel to produce y
		# Generate random channel noise and then received signal y
		z = np.random.randn(n, 1) * sigma
		y = (x + z).reshape(-1, 1)

		# use an initialisation with 50-70% of sections decoded correctly and the other sections set to zero
		beta_0 = beta/np.sqrt(n*P/L)
		beta_0[:L_zero*M] = 0


		# hard decision initialisation
		# peel off contributions from the decoded sections from y to generate y'
		x_remove = Ab(beta_0)
		y_new = y - x_remove
		# run amp on the remaining sections
		Ab_new, Az_new = sparc_transforms_shorter(L_zero, M, n, ordering)

		beta_hard = amp(y_new, sigma, Pl[:L_zero], L_zero, M, T, Ab_new, Az_new).reshape(-1)
		# Convert decoded beta_hard back to a message
		# beta_hard gives you the first L_zero sections of the final received message. 
		rx_message_hard = []
		for l in range(L_zero):
			idx = np.argmax(beta_hard[l*M:(l+1)*M])
			rx_message_hard.append(idx)
		# find the errors in this rx message (note that we don't need to check for errors in 
		# the sections we said were decoded correctly and we know there won't be any.)
		# compute BER for this hard initialisation method
		ber_hard = ber_hard + sum(bin(a^b).count('1') for (a, b) in zip(sparc_indices[:L_zero], rx_message_hard))/total_bits      
		
		

		# Soft initialisation
		# initialise with beta_0 and the original channel output y.
		# look at the error performance of this. 
		beta_soft = amp(y, sigma, Pl, L, M, T, Ab, Az, beta_0).reshape(-1)
		rx_message_soft = []
		for l in range(L):
			idx = np.argmax(beta_soft[l*M:(l+1)*M])
			rx_message_soft.append(idx)
		# compute BER for the soft intialisation method
		ber_soft = ber_soft + sum(bin(a^b).count('1') for (a, b) in zip(sparc_indices, rx_message_soft))/total_bits 

		# compare both to amp decoding with no initialisation
		beta_no_init = amp(y, sigma, Pl, L, M, T, Ab, Az).reshape(-1)   
		rx_message_no_init = []
		for l in range(L):
			idx = np.argmax(beta_no_init[l*M:(l+1)*M])
			rx_message_no_init.append(idx)
		# compute BER for no initialisation
		ber_no_init = ber_no_init + sum(bin(a^b).count('1') for (a, b) in zip(sparc_indices, rx_message_no_init))/total_bits
	ber_hard = ber_hard/repeats
	ber_soft = ber_soft/repeats
	ber_no_init = ber_no_init/repeats

	print("The ber for hard init is: ", ber_hard)
	print("The ber soft init is: ", ber_soft)
	print("The ber for no init is: ", ber_no_init)

