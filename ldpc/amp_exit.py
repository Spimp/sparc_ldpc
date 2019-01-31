# Simulating the effect of varying I_A on the output of the AMP decoder to produce an EXIT chart
import numpy as np
import math as math
from pylab import * 
import py.ldpc as ldpc
import matplotlib.pyplot as plt
import time
import csv
from bitarray import bitarray
from sparc_ldpc import * 

class imported_E:
	def __init__(self, X, E, I_a, SNR_dB):
		self.X = X
		self.E = E
		self.I_a = I_a 
		self.SNR_dB = SNR_dB

# J and J_inverse based on approximation in the appendix of Design of LDPC Codes for modulation and detection.
def J_inverse(I):
	assert I>=0 and I<=1
	if (I<=0.3646):
		return 1.09542*(I**2) + 0.214217*(I) + 2.33727*np.sqrt(I)
	else:
		return -0.706692*log(0.386013*(1-I)) + 1.75017*(I) 

def J(sigma):
	assert sigma>=0
	if (sigma<=1.6363):
		return -0.0421061*(sigma**3) + 0.209252*(sigma**2) + -0.00640081*(sigma)
	elif (sigma<10):
		return 1- np.exp(0.00181491*(sigma**3) - 0.142675*(sigma**2) - 0.0822054*sigma + 0.0549608)
	else:
		return 1

# generate a random array populated with -1s and 1s 
def gen_bits(length):
	# multiply by -2 and add one so have -1s and 1s
	return (np.random.randint(0,2,length)*-2)+1

# Function produces y from X a series of bits for inputting into the amp() function
# X: array of LlogM input bits
# L, M, N are all SPARC parameters
# sigma_w is the variance of the noise
# P is the power that will be uniformly allocated
# returns: y, Ab, Az, and Pl
def prep_y(X, L, M, n, sigma_w, P):
	# Calculate the power allocation
	Pl = P/L * np.ones(L)
	# Need to convert X, which contains -1 and +1s to bits
	X = (X-1)*-1/2
    # convert the bits to indices for encoding
	X_indices = bits2indices(X, M)

	assert len(X_indices)==L
	   
	# Generate the SPARC transform functions A.beta and A'.z
	Ab, Az, ordering = sparc_transforms(L, M, n)

	# Generate our transmitted signal x
	β_0 = np.zeros((L*M, 1))
	for l in range(L):
	    β_0[l*M + X_indices[l]] = np.sqrt(n * Pl[l])
	x = Ab(β_0)

	# Generate random channel noise and then received signal y
	w = np.random.randn(n, 1) * sigma_w
	y = (x + w).reshape(-1, 1)

	return y, Ab, Az, Pl

# a and b are 2 numpy array of the same length.
# When there is a zero in both a and b at the same index, this is removed from both arrays. 
def remove_common_zeros(a, b):
	# find the indices of all zeros in a
	a_indices = np.where(a==0)[0]
	# find the indices of all zeros in b
	b_indices = np.where(b==0)[0]
	# Find all the indices that appear in both
	remove_indices = np.intersect1d(a_indices, b_indices)
	# remove these indices from both arrays and return
	a = np.delete(a, remove_indices, None)
	b = np.delete(b, remove_indices, None)
	return a, b


# calculate E using X, I_a, snr_dB and the SPARC parameters
# X is the randomly generated bits
# I_a is the mutual information of the a priori input to the AMP decoder
# snr_dB is the snr of the simulation of the channel 
# E is the LLRs on the a posteriori information from the AMP decoder
# returns E
def calc_E(X, I_a, snr_dB, SPARCParams, csv_filename=None):
	# Sparc parameters
	L = sparcparams.L
	M = sparcparams.M
	logm = np.log2(M)
	P=sparcparams.p
	r_sparc = sparcparams.r
	T = sparcparams.t
	n = int(L*np.log2(M) / r_sparc)
	snr = 10**(snr_dB/20)
	sigma_w = np.sqrt(P/snr)

	PE_pos = 0
	PE_neg = 0
	bin_edges_pos = 0
	bin_edges_neg = 0
	#calculate sigma_a and mu_a
	sigma_a = J_inverse(I_a)
	mu_a = (sigma_a**2)/2
	#print("mu_a: ", mu_a)
	# N_a is gaussian noise with zero mean and variance sigma_a**2. We want an array of samples from this the same size as X. 
	N_a = np.random.randn(len(X))*sigma_a
	#print("N_a: ", N_a)
	# the simulated L-values coming out of the LDPC decoder
	A = mu_a*X + N_a 
	# convert LLRs in A to bitwise posterior probabilities
	bitwise_a = 1/(1+np.exp(A))
	# convert this to sectionwise posterior probabilities
	beta_0 = bp2sp(bitwise_a, L, M)
	# perform amp decoding on these beta_0
	sigma=None # this value isn't actually used in the function so doesn't matter
	# Generate y, Ab, and Az to pass into amp()
	y, Ab, Az, Pl = prep_y(X, L, M, n, sigma_w, P)
	beta_T = amp(y, sigma, Pl, L, M, T, Ab, Az, beta_0*np.sqrt(n*P/L)).reshape(-1)

	# convert sectionwise posterior probabilities to bitwise posterior probabilities
	sectionwise = beta_T/np.sqrt(n*P/L)
	bitwise_e = sp2bp(sectionwise, L, M)
	# clip the bitwise_e to avoid a divide by zero error in the log
	#np.clip(bitwise_e, 0.00000000000000000000000001, 1-0.00000000000000000000000001, out=bitwise_e)
	# convert bitwise post. to LLRs 
	E = np.log(1-bitwise_e)-np.log(bitwise_e)
	# set -inf and +inf to real numbers with v large magnitude
	np.nan_to_num(E)
	#print("E: ", E)

	if csv_filename!=None:
		myFile = open(csv_filename, 'a')
		with myFile:
			myFields = ['I_a', 'snr_dB','X', 'E']
			writer = csv.DictWriter(myFile, fieldnames=myFields)
			writer.writeheader()
			writer.writerow({'I_a':I_a, 'snr_dB':snr_dB,'X':X,'E':E})

	return E

# given E and X return the histograms of E
# X = generated input bits
# bin_number, max_bin, min_bin determine the number of histogram bins and the range of the bins
# plot is a bool. Set true to plot histograms
# Return: histograms of E, and their means and variances. 
def hist_E(X, E, bin_number=500, max_bin=40, min_bin=-40, plot=False):
	assert(len(E)==len(X))

	# find the indices of all the occurences of a 1 in X
	index_pos1 = np.where(X==1)[0]

	# find the indices of all the occurences of a -1 in X
	index_neg1 = np.where(X==-1)[0]

	# calculate histogram of P(E(X_i)|X_i=+1). 
	bin_width = (max_bin-min_bin)/(bin_number-1)
	# set the bin edges so that both PEs are over the same range and have the same width bins. 
	bin_edges = np.linspace(min_bin, max_bin, bin_number)
	PE_pos, bin_edges_pos = np.histogram(E[index_pos1], bins=bin_edges, density=True)

	PE_neg, bin_edges_neg = np.histogram(E[index_neg1], bins=bin_edges, density=True)

	# the midpoints of the bins
	mids = 0.5*(bin_edges[1:]+bin_edges[:-1])
	# the means for the two histograms. Calculated by doing a weighted average of the mids
	mean_pos = np.average(mids, weights=PE_pos)
	mean_neg = np.average(mids, weights=PE_neg)
	# the variance for the two histograms
	var_pos = np.average((mids - mean_pos)**2, weights=PE_pos)
	var_neg = np.average((mids - mean_neg)**2, weights=PE_neg)
	if plot:
		# Plot histogram of P(E(X_i)|X_i=+1).
		plt.figure(1)
		plt.hist(E[index_pos1], bins=bin_edges, density = True)
		plt.title("Normalised histogram of $P(E(X_i)|X_i=+1)$. $\mu_a=$"+str(mu_a)+" $SNR_{dB}=$"+str(snr_dB))

		plt.figure(2)
		plt.hist(E[index_neg1], bins=bin_edges, density = True)
		plt.title("Normalised histogram of $P(E(X_i)|X_i=-1)$. $\mu_a=$"+str(mu_a)+" $SNR_{dB}=$"+str(snr_dB))
		plt.show()

	return PE_pos, PE_neg, mean_pos, mean_neg, var_pos, var_neg, bin_width

# calculate I_e from the relevant histograms of E and the bin_width used for these histograms
def calc_I_e(PE_pos, PE_neg, bin_width):
	# remove all the common zeros to avoid divide by zero. 
	PE_pos, PE_neg = remove_common_zeros(PE_pos, PE_neg)

	# produce what is essentially a new histogram for each E that corresponds to the contents of the integral
	integral_neg = PE_neg * log2(2*PE_neg/(PE_neg+PE_pos))
	integral_pos = PE_pos * log2(2*PE_pos/(PE_neg+PE_pos))
	
	# when PE_neg = 0 or PE_pos=0 we end up with 0*log2(0) which in python evaluates to nan
	# Below we replace all the NaNs with a zero. This could either be 0 or -1 depending on our convention.
	integral_neg[np.isnan(integral_neg)] = 0
	integral_pos[np.isnan(integral_pos)] = 0
	# calculate I_e by summing over the contents of the integral and summing for each value of x to approximate the equation for I_e
	I_e = 1/2 * (bin_width*sum(integral_neg) + bin_width*sum(integral_pos))
	#print(I_e)
	return I_e

# import E and the corresponding X, I_a, and snr_dB from file.
# store them in a new class which has a matrix of the repeats of E
# and of X and the corresponding I_a and snr_dB
# return a dictionary of these object, each indexed by str(I_a)+' '+str(SNR_dB) 
def import_E_fromfile(fileName, datapoints, repeats, Llogm):
	imported_E_dict = {}
	with open(fileName) as myfile:		
		for j in range(10):
			E = np.zeros((repeats, Llogm))
			X = np.zeros((repeats, Llogm))
			SNR_dB = np.zeros(repeats)
			I_a = np.zeros(repeats)
			i=0
			for i in range(repeats):
				reader = csv.DictReader(myfile)
				E[i,:] = row['E']
				X[i,:] = row['X']
				I_a[i] = row['I_a']
				SNR_dB[i] = row['snr_dB']
			imported_E_dict[str(I_a[0])+' '+str(SNR_dB[0])] = imported_E(X=X, E=E, I_a=I_a[0], SNR_dB=SNR_dB[0])
	return imported_E_dict


if __name__ == "__main__":
	# NOTE have set T=10!!!!
	t0=time.time()
	L=768
	M=512
	logm = np.log2(M)
	sparcparams = SPARCParams(L=L, M=M, sigma=None, p=1.8, r=0.877, t=10)

	'''
	# just plotting one set of histograms
	X = gen_bits(int(L*logm))
	print(X)
	snr_dB = 8.0
	I_a = 0.4

	E = calc_E(X, I_a, snr_dB, sparcparams)
	PE_pos, PE_neg, mean_pos, mean_neg, var_pos, var_neg, _ = hist_E(X, E, bin_number=500, max_bin=40, min_bin=-40, plot=True)
	print("positive mean: ", mean_pos)
	print("negative mean: ", mean_neg)
	print("positive variance: ", var_pos)
	print("negative variance: ", var_neg)
	'''
	# plotting the EXIT chart for the AMP decoder for a range of SNR
	repeats = 10
	datapoints = 5
	I_a_range = np.linspace(0, 0.9, 10)
	P = 1.8
	for k in range(repeats):
		snr_dB = np.linspace(6, 10, datapoints)
		# accumulative values of I_e for each snr value
		I_e_accum = np.zeros((datapoints,10))
		j=0
		for s_dB in snr_dB:
			I_e = np.zeros(10)
			i=0
			for I_a in I_a_range:
				# randomly generate bits
				X = gen_bits(int(L*logm))
				#print(X)

				# generate the histograms for E and some statistics about them
				E = calc_E(X, I_a, s_dB, sparcparams, csv_filename='E_data_L768_M512_r0_877_p1_8_t10.csv')

				PE_pos, PE_neg, mean_pos, mean_neg, var_pos, var_neg, bin_width = hist_E(X, E, bin_number=500, max_bin=40, min_bin=-40, plot=False)
	
				# calculate I_e
				I_e[i] = calc_I_e(PE_pos, PE_neg, bin_width)
				i=i+1
			I_e_accum[j,:] = I_e_accum[j,:] + I_e
			j=j+1
	I_e_accum = I_e_accum/repeats
	print("Wall clock time elapsed: ", time.time()-t0)

	fig, ax = plt.subplots()
	ax.plot(I_a_range, I_e_accum[0,:], 'b--', label='$SNR_{dB}$='+str(snr_dB[0]))	
	ax.plot(I_a_range, I_e_accum[1,:], 'k--', label='$SNR_{dB}$='+str(snr_dB[1]))
	ax.plot(I_a_range, I_e_accum[2,:], 'm--', label='$SNR_{dB}$='+str(snr_dB[2]))
	ax.plot(I_a_range, I_e_accum[3,:], 'c--', label='$SNR_{dB}$='+str(snr_dB[3]))
	ax.plot(I_a_range, I_e_accum[4,:], 'r--', label='$SNR_{dB}$='+str(snr_dB[4]))
	
	plt.xlabel('$I_A$')
	plt.ylabel('$I_E$')
	plt.legend(loc=6, prop={'size': 7})
	plt.title("The EXIT chart for the AMP decoder")
	plt.savefig('amp_exitchart_L768_M512_10reps_t10.png')	
	
	