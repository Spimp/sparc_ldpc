# Simulating the effect of varying I_A on the output of the AMP decoder to produce an EXIT chart
import numpy as np
import math as math
from pylab import * 
import py.ldpc as ldpc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
import csv
from bitarray import bitarray
from sparc_ldpc import * 
import EXIT_chart as exit

class imported_E:
	"""
	class used to store E when it is imported from a csv file. 
	"""
	def __init__(self, X, E, I_a, SNR_dB):
		self.X = X
		self.E = E
		self.I_a = I_a 
		self.SNR_dB = SNR_dB
	def __eq__(self, other):
		# check if two sets of E data are the same. 
		return (self.X==other.X).all() and (self.E==other.E).all() and self.I_a==other.I_a and self.SNR_dB==other.SNR_dB

# J and J_inverse based on approximation in the appendix of Design of LDPC Codes for modulation and detection.
def J_inverse(I):
	assert I>=0 and I<=1
	if I==1:
		I=np.clip(I,a_min=None, a_max=0.9999)
		print("Warning clipping I from 1 to 0.9999")
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
 
# This function takes in our sectionwise posterior probabilities as beta, which is the output
# from the bp2sp conversion performed on the LLRS from A. 
# it then produces a hard initialisation for the amp decoder y_new which only contains sections
# which were not decoded well in the A and then bp2sp conversion.
def hard_initialisation(beta, L, M, n, ordering, y, Pl, Ab, threshold=0.5, ldpc_sections=None):
	"""
	Used in threshold initialised information exchange to produce the reduced channel output and the 
	functions for the smaller design matrix which only contains sections which didn't exceed the threshold. 
	beta: sectionwise posterior probabilities output from the bp2sp conversion. 
	L: number of sections in beta
	M: number of columns per section
	ordering: the ordering of the design matrix A. Allows you to select the relevant columns of the design matrix
	y: the original channel output
	Pl: Vector containing the power allocation. Should have one entry per section of beta
	Ab: the function used to apply the design matrix to a beta vector
	threshold: above this outputs from LDPC decoder are hard decoded
	ldpc_sections: number of sections covered by the LDPC code. 
	Returns:
	y_new: reduced channel output
	Ab_new & Az_new: reduced functions used for the design matrix. Just for sections that don't exceed threshold
	amp_sections: sections which should have amp decoding performed on them
	L_amp_sections: the number of sections in amp_sections.
	"""
	if ldpc_sections==None:
		# all sections covered by the ldpc code
		ldpc_sections = L
	# note that beta is already normalized
	beta_0 = beta

	# store the sections which we still want to perform amp decoding on
	amp_sections = []
	# keep track of the number of sections in the amp
	L_amp_sections = 0
	for l in range(L):
		if l<L-ldpc_sections:
			# if the section is not covered by the LDPC code, we always want it to be involved in the amp decoding. 
			idx = np.array([[],[]])
		else:
			# find the index of any element that has probability greater than threshold
			idx = np.where(beta_0[l*M:(l+1)*M]>threshold)
			#assert idx[0].size<=1
		# if the section contains a single position with the probability greater than the threshold
		if idx[0].size==1:
			# set all positions to zero
			beta_0[l*M:(l+1)*M] = 0
			# set the position corresponding to the probability greater than the threshold to the power allocation 
			# NOTE: previously I was setting this to 1 which was an error. 
			beta_0[(l*M)+idx[0]] = np.sqrt(n * Pl[l])
		else:
			# else the section does not contain any high probability components
			# set all sections to zero
			beta_0[l*M:(l+1)*M]=0
			# add the section to amp_sections so we can include it in the rows kept from ordering
			amp_sections.append(l)
			L_amp_sections += 1
		idx=None
	# hard decision initialisation
	# peel off contributions from the decoded sections from y to generate y'

	x_remove = Ab(beta_0)
	y_new = y - x_remove
	ordering_reduced = ordering[amp_sections,:]
	if L_amp_sections>0:
		# run amp on the remaining sections
		Ab_new, Az_new = sparc_transforms_shorter(L_amp_sections, M, n, ordering_reduced)
	else:
		Ab_new = None
		Az_new = None
	#print("L_amp_sections: ", L_amp_sections)

	return y_new, Ab_new, Az_new, amp_sections, L_amp_sections


def prep_y(X, L, M, n, sigma_w, P, a=None, f=None, C=None):
	"""
	Produces the channel output y from the input X 
	X: array of LlogM input bits taking values -1 and +1 
	L, M, N are all SPARC parameters
	sigma_w: is the standard deviation of the noise
	P: total power 
	returns: y, Ab, Az, and Pl
	"""
	# Calculate the power allocation
	if a==None:
		# uniform power allocation
		Pl = P/L * np.ones(L)
	else:
		# parameterised power allocation (fn from sprac_ldpc.py)
		Pl = pa_parameterised(L, C, P, a, f)
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
	#print("If this close to zero, power allocation is correct: ", (sum(x**2)/n)-P)

	# Generate random channel noise and then received signal y
	w = np.random.randn(n, 1) * sigma_w
	y = (x + w).reshape(-1, 1)

	return y, Ab, Az, Pl, ordering

def remove_common_zeros(a, b):
	"""
	a and b are 2 numpy array of the same length.
	When there is a zero in both a and b at the same index, this is removed from both arrays.
	"""
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
def calc_E(X, I_a, snr_dB, sparcparams, csv_filename=None, threshold=0.5):
	"""
	calculates E
	X: randomly generated bits taking values -1 and 1
	I_a: the mutual information of the a priori input to the AMP decoder
	snr_dB: the snr of the simulation of the channel
	Returns E: the extrinsic log likelihood ratios (LLRs) produced by the AMP decoder
	"""
	# Sparc parameters
	L = sparcparams.L
	M = sparcparams.M
	logm = int(np.log2(M))
	P=sparcparams.p
	r_sparc = sparcparams.r
	T = sparcparams.t
	a = sparcparams.a
	f = sparcparams.f
	C = sparcparams.C

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
	#print("max in beta_0: ", np.max(beta_0))
	# perform amp decoding on these beta_0
	sigma=None # this value isn't actually used in the function so doesn't matter
	# Generate y, Ab, Pl, and ordering to pass into hard)initialisation()
	y, Ab, _, Pl, ordering = prep_y(X, L, M, n, sigma_w, P, a, f, C)

	y_new, Ab_new, Az_new, amp_sections, L_amp_sections = hard_initialisation(beta_0, L, M, n, ordering, y, Pl, Ab, threshold)
	#print("L_amp_sections: ", L_amp_sections)
	# keep the LLRs corresponding to the sections hard decoded before performing amp
	E = A
	# if there are no amp sections we just keep E=A.
	if L_amp_sections>0:
		beta_T = amp(y_new, sigma, Pl[amp_sections], L_amp_sections, M, T, Ab_new, Az_new).reshape(-1)

		# convert sectionwise posterior probabilities to bitwise posterior probabilities
		sectionwise = beta_T/np.sqrt(n*np.repeat(Pl[amp_sections],M))
		#print("sum of sectionwise is: ", np.sum(sectionwise))
		#print("sum of sectionwise should be: ", L_amp_sections)
		bitwise_e = sp2bp(sectionwise, L_amp_sections, M)
		# convert bitwise post. to LLRs 
		E_amp_sections = np.log(1-bitwise_e)-np.log(bitwise_e)
		# set -inf and +inf to real numbers with v large magnitude
		E_amp_sections=np.nan_to_num(E_amp_sections)
	
		# convert amp_sections to a list of all the positions in these sections
		amp_positions = np.tile(np.linspace(0,logm-1,logm,dtype=np.dtype(np.int16)),L_amp_sections)
		amp_positions = amp_positions + logm*np.repeat(amp_sections, logm)
		# replace the LLRs for the amp_sections with the new LLRs
		E[amp_positions] = E_amp_sections
	###########################
	#Note testing out line below. CHECKED and seems to work fine. 
	# doing this so that values that get set to inf get counted in the histogram 
	# But actually I don't think this was happening much so this line is maybe redundant. 
	np.clip(E, -55, 55, out=E)
	#print("E: ", E)
	np.set_printoptions(threshold=np.nan)
	if csv_filename!=None:
		myFile = open(csv_filename, 'a')
		with myFile:
			myFields = ['I_a', 'snr_dB','X', 'E']
			writer = csv.DictWriter(myFile, fieldnames=myFields)
			writer.writeheader()
			writer.writerow({'I_a':I_a, 'snr_dB':snr_dB,'X':X,'E':E})

	return E

def hist_E(X, E, bin_number=500, max_bin=40, min_bin=-40, plot=False, snr_dB='Not given'):
	"""
	Produce histogram of E the extrinsic LLRs
	X: generated input bits 
	E: extrinisic LLRs 
	bin_number: number of bins used for the histogram
	max_bin, min_bin: set the range over which the histogram is plotted
	plot: bool to determine if this function plots the histogram values it is producing
	snr_dB: the snr in dB that was used to produce E. This is only used for the figure headings 
	Return: 
	PE_pos: the histogram/ probability distribution of E values which have positive input X
	PE_neg: the histogram of E values for negative input X
	mean_pos: the mean of the positive histogram
	mean_neg: the mean of the negative histogram
	var_pos: the variance of the positive histogram
	var_neg: the variance of the negative histogram
	bin_width: the width of each bin
	"""
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
		plt.title("Normalised histogram of $P(E(X_i)|X_i=+1)$."+" $SNR=$"+str(snr_dB)+"dB")

		plt.figure(2)
		plt.hist(E[index_neg1], bins=bin_edges, density = True)
		plt.title("Normalised histogram of $P(E(X_i)|X_i=-1)$."+" $SNR=$"+str(snr_dB)+"dB")
		plt.show()

	return PE_pos, PE_neg, mean_pos, mean_neg, var_pos, var_neg, bin_width

def calc_I_e(PE_pos, PE_neg, bin_width):
	"""
	Calculate I_e, the extrinsic mutual information 
	Implementing the integral shown in step g section 4.3.1 of the report
	PE_pos: histogram of E values for positive input X
	PE_neg: histogram of E values for negative input X
	bin_width: the width of the histogram bins
	return: I_es
	"""
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
 
def import_E_fromfile(fileName, datapoints, repeats, Llogm):
	"""
	If data from a previous run has been exported to a csv file, this function reads that data in 
	imports E and the corresponding X, I_a, and snr_dB from file.
	store them in a new class which has a matrix of the repeats of E
	and of X and the corresponding I_a and snr_dB
	return a dictionary of these object, each indexed by str(I_a)+' '+str(SNR_dB)
	fileName: a csv filename
	datapoints: the number of snrs included in the data
	repeats: the number of times each datapoint is repeated
	Llogm: the number of sections multiplied by the log of the number of columns in each section (i.e. the total number of bits represented by the SPARC)
	"""

	# allows you pick out numbers from strings. 
	numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
	rx = re.compile(numeric_const_pattern, re.VERBOSE)
	imported_E_dict = {}
	with open(fileName) as myfile:
		for i in range(repeats):
			for k in range(10):
				for j in range(datapoints):
					E = np.zeros(Llogm)
					X = np.zeros(Llogm)
					SNR_dB = 0
					I_a = 0
					reader = csv.DictReader(myfile)
					for row in reader:	# this should only actually loop through one row before it gets stuck?
						E = np.float_(rx.findall(row['E']))
						X = np.int_(rx.findall(row['X']))
						I_a = float(row['I_a'])
						SNR_dB = float(row['snr_dB'])
						break
					# this ensures that each time a result is repeated, it is stored under a new key.
					# Results don't always appear to be output in order so you can't always just use i for this. 
					for l in range(repeats):
						if str(np.round(I_a,1))+' '+str(int(np.round(SNR_dB)))+' '+str(int(l)) in imported_E_dict:
							if imported_E_dict[str(np.round(I_a,1))+' '+str(int(np.round(SNR_dB)))+' '+str(int(l))] == imported_E(X=X, E=E, I_a=I_a, SNR_dB=SNR_dB):
								print("repeat dict")
								break
							else:
								continue
						else:
							imported_E_dict[str(np.round(I_a,1))+' '+str(int(np.round(SNR_dB)))+' '+str(int(l))] = imported_E(X=X, E=E, I_a=I_a, SNR_dB=SNR_dB)
							break
					#print(str(np.round(I_a,1))+' '+str(int(np.round(SNR_dB)))+' '+str(int(i)))
	return imported_E_dict

def polynomial(I_a, I_e):
	"""
	Produces a 3rd order polynomial to fit the exit curve
	I_a: array of apriori mutual information
	I_e: array of the corresponding extrinsic mutual information
	return: c array of polynomial coefficients
	c[0]: the constant
	c[3]: is the coefficient of I_a^3
	"""
	a = np.zeros((len(I_a), 4))
	for i in range(4):
		a[:,i] = I_a**i
	c = np.linalg.lstsq(a, I_e, rcond=-1)[0]
	return c	


# See section 4.4 Combined EXIT chart in the report
def I_E_VND_amp(I_A_VND, d_v, poly_coeff):
	"""
	Calculates the value of the extrinsic mutual information for the combined AMP and VND decoder for the input I_A_VND
	I_A_VND: a priori mutual information for the VND decoder. Single value, does not take an array
	d_v: degree of the variable node in question
	poly_coeff: the coefficients of the 3rd order polynomial used to approx the amp exit curve
	return: the extrinsic mutual information output from the VND. 
	"""
	# calculate a priori mutual information for the AMP decoder
	I_A_amp = J(np.sqrt(d_v)*J_inverse(I_A_VND))
	# calculate the extrinsic mutual information for the AMP decoder using the 3rd order polynomial approximation
	I_E_amp = np.sum(poly_coeff * np.array([1, I_A_amp, I_A_amp**2, I_A_amp**3]))
	# clip I_E_amp at a max of I_E_amp=0.9999 as it must be less than 1
	I_E_amp = np.clip(I_E_amp, a_min=None, a_max=0.9999)
	I_E_VND = J(np.sqrt((d_v-1)*(J_inverse(I_A_VND))**2+(J_inverse(I_E_amp))**2))
	return I_E_VND

def I_E_VND_amp_array(I_A_VND_array, a_v, b_v, poly_coeff):
	"""
	Calculate the array of extrinsic mutual information for the combined AMP and VND decoder. 
	I_A_VND_array: array of a priori mutual information input to the VND decoder
	a_v: fraction of variable nodes with degree equal to its indes
	b_v: fraction of edges incident to variables nodes of degree equivalent to the index 
	poly_coeff: the coefficients of the 3rd order polynomial used to approx the amp exit curve
	"""
	I_E_VND = np.zeros(I_A_VND_array.shape)
	for i in range(len(a_v)):
		# a_v gives the fraction of variable nodes with degree equal to its index.
		# when this is zero means there are no v nodes with this degree so skip it
		if a_v[i]!=0:
			for j in range(len(I_E_VND)):
				I_E_VND[j] = I_E_VND[j]+b_v[i]*I_E_VND_amp(I_A_VND_array[j], i, poly_coeff)
	return I_E_VND


def plot_amp_ldpc_exit(poly_coeff, a_v, d_c, R_ldpc):
	"""
	Plot the exit chart for the combined SPARC and LDPC code
	AMP decoder and VND are combined into one curve
	CND is the other curve
	Allows you to see if the two curves intersect
	This function is for a regular check node decoder 
	poly_coeff: coefficients of the 3rd ordder polynomial representing the amp exit curve 
	a_v: array giving the fraction of variable nodes of the degree equal to its index
	d_c: the degree of the check node. Should just be one value.
	R_ldpc: the rate of the LDPC code 
	"""
	dv_degrees = np.linspace(0,len(a_v)-1,len(a_v),dtype=np.dtype(np.int16))
	# calculating the fraction of edges incident to variables nodes of degree equivalent to the index
	b_v = (dv_degrees*a_v)/((1-R_ldpc)*d_c)

	# plot I_E,VND as a func of I_A,VND and I_A_amp
	I_A_VND = linspace(0,1,21)
	I_E_VND = I_E_VND_amp_array(I_A_VND, a_v, b_v, poly_coeff)

	I_A_CND = linspace(0.1, 1, 19)
	I_E_CND = np.zeros(I_A_CND.shape)
	for k in range(len(I_E_CND)):
		I_E_CND[k] = 1 - exit.I_E_REP(I_A_CND[k], d_c)
	fig, ax = plt.subplots()
	plt.plot(I_A_VND, I_E_VND, color='k', linestyle='-', label = '$VND_{comb}$')
	plt.plot(I_E_CND, I_A_CND, color='r', linestyle='-', label = '$CND$')
	plt.xlabel('$I_{A,VND_{comb}}, I_{E,CND}$', fontsize=18) 
	plt.ylabel('$I_{E,VND_{comb}}, I_{A,CND}$', fontsize=18)
	plt.tight_layout()
	plt.legend()
	plt.show()

# This function, along with trial and error was used to find valid degree distributions
def plane_intersect(a, b, plot=False):
	"""
	a, b   4-tuples/lists
	       Ax + By +Cz + D = 0
	       A,B,C,D in order  

	output: a point on line of intersection and the direction of the line of intersection
	"""
	a_vec, b_vec = np.array(a[:3]), np.array(b[:3])
	# direction of line of intersection
	aXb_vec = np.cross(a_vec, b_vec)
	# point on line of intersection
	A = np.array([a_vec[:2],b_vec[:2]])
	x = np.linalg.solve(A, np.array([-1*a[3],-1*b[3]]))
	# so point on plane is [x[0],x[1],0]
	point = np.array([x[0],x[1],0])
	if plot:	
		# find a1 and a2 in terms of a3
		a3 = np.array([0,1])
		t = aXb_vec[2]*a3
		a1 = x[0] + aXb_vec[0]*t
		a2 = x[1] + aXb_vec[1]*t
		print("A point on the line is:", point)
		print("The direction of the line is: ", aXb_vec)
		fig, ax = plt.subplots()
		plt.plot(a1,a2,'k-')
		plt.axhline(y=0, color='r', linestyle='-', label='Shannon limit')
		plt.plot()
		plt.xlabel("$a_1$")
		plt.ylabel("$a_2$")
		plt.show()
	return point, aXb_vec


def amp_exit_curve(sparcparams, low_snr_dB, high_snr_dB, repeats, x_axis_points, threshold, poly_curve=0, bin_number=500, import_data=False, export_csv_filename=None, import_csv_filename=None):
	"""
	This functions plots the AMP decoder exit curves for 4 different SNRs. 
	This function only works for threshold initialised information exchange. 
	It then finds a 3rd order polynomial fit for one of these curves and prints this
	It plots this polynomial curve with the curve it is trying to represent
	low_snr_dB: the lowest snr in decibels that an exit curve is plotted for 
	high_snr_dB: the highest snr in decibels that an exit curve is plotted for 
	import_data: bool to determine if the data should be imported. This speeds up this function
	export_csv_filename: name of the csv file that data is exported to. If this is None then no data will be exported to csv
	import_csv_filename: name of the csv file that data is imported from 
	repeats: the number of times each datapoint is repeated
	x_axis_points: the number of points on the x axis. I.e. the number of different I_as that are plotted 
	sparcparams: the parameters of the sparc code being used
	bin_number: the number of bins used for plotting the histograms of E
	poly_curve: the curve which is used to generate the polynomial fit. This must be between 0 and curves-1
	threshold: the threshold used in threshold initialised information exchange. 
	"""
	# SPARC parameters
	L = sparcparams.L
	M = sparcparams.M
	logm = np.log2(M)
	P=sparcparams.p
	R_sparc = sparcparams.r
	t = sparcparams.t
	# the number of curves that are plotted/ the number of snrs included
	# 4 curves will always be plotted unless this is changed. 
	# Note: if want to change this, also need to change code for plotting figure 1
	curves=4
	I_a_range = np.linspace(0, 0.99, x_axis_points)
	if import_data==True:	# if not exporting, want to import E into a dictionary
		if import_csv_filename==None:
			print("Please enter a valid csv filename in import_csv_filename")
		# sometimes need to set the repeats in here higher than the actual repeats to ensure all the data is imported. I don't understand why!
		imported_E_dict	= import_E_fromfile(fileName = import_csv_filename, datapoints = curves, repeats = repeats, Llogm = int(L*logm))
		print(len(imported_E_dict))
	# accumulative values of I_e for each snr value
	I_e_accum = np.zeros((curves,x_axis_points))
	# work in snr for the EXIT charts as then don't have to work about EbN0 and rate. 
	snr_dB = np.linspace(low_snr_dB, high_snr_dB, curves)
	for k in range(repeats):
		j=0
		for s_dB in snr_dB:
			I_e = np.zeros(x_axis_points)
			i=0
			for I_a in I_a_range:
				a = None
				if import_data==False: # if import_data is true, we already have these values.
					# randomly generate bits
					X = gen_bits(int(L*logm))
					#print(X)

					# generate the histograms for E and some statistics about them
					# if export_csv_filename is None, data won't be exported to a csv file.
					E = calc_E(X, I_a, s_dB, sparcparams, threshold=threshold, csv_filename=export_csv_filename)
				else:	
					# get the required entry by using a key which is 'I_a s_dB k' where k is the current repetition
					a = imported_E_dict[str(np.round(I_a,1))+' '+str(int(np.round(s_dB)))+' '+str(k)]

					X = a.X
					assert(len(X)==int(L*logm))
					E = a.E
					if(len(E)!=int(L*logm)): print("i: ", i, "I_A:", I_a, "j: ", j)
					assert(len(E)==int(L*logm))
				#print('X[0] is: ', X[0])
				#print("E[0] is: ", E[0])
				# produce histograms of E
				PE_pos, PE_neg, mean_pos, mean_neg, var_pos, var_neg, bin_width = hist_E(X, E, bin_number=bin_number, max_bin=60, min_bin=-60)#, plot=True, snr_dB=s_dB)
	
				# calculate I_e using the histograms of E 
				I_e[i] = calc_I_e(PE_pos, PE_neg, bin_width)
				i=i+1
			I_e_accum[j,:] = I_e_accum[j,:] + I_e
			j=j+1
	# Average I_e over the repeats
	I_e_accum = I_e_accum/repeats
	# calculate the polynomial coefficients for the 3rd order polynomial representation of the exit curve
	poly_coeff = polynomial(I_a_range, I_e_accum[poly_curve,:]) 
	# Note that this will give the constant first and the coefficient of I_a**3 last
	print("The coefficients for the polynomial are: ",poly_coeff)
	print("Wall clock time elapsed: ", time.time()-t0)

	# plot the exit curves
	fig, ax = plt.subplots()
	ax.plot(I_a_range, I_e_accum[0,:], 'b--', label='$SNR$='+str(snr_dB[0])+'$dB$')	
	ax.plot(I_a_range, I_e_accum[1,:], 'k--', label='$SNR$='+str(snr_dB[1])+'$dB$')
	ax.plot(I_a_range, I_e_accum[2,:], 'm--', label='$SNR$='+str(snr_dB[2])+'$dB$')
	ax.plot(I_a_range, I_e_accum[3,:], 'c--', label='$SNR$='+str(snr_dB[3])+'$dB$')
	#ax.plot(I_a_range, I_e_accum[4,:], 'r--', label='$SNR$='+str(snr_dB[4])+'$dB$')
	plt.xlabel('$I_A$')
	plt.ylabel('$I_E$')
	plt.legend(loc=6, prop={'size': 7})
	#plt.title("The EXIT chart for the AMP decoder")
	#plt.savefig('amp_exitchart_L128_M4_40reps_500bins_r1_5_P2.png')	
	#plt.savefig('amp_exit_threshold_L256M32R1P4Bins350Threshold0_99_200reps.pdf')	
	
	# plot the polynomial curve fitted to the relevant exit curve
	fig, ax = plt.subplots()
	ax.plot(I_a_range, I_e_accum[poly_curve,:], 'k--', label='$SNR$='+str(snr_dB[poly_curve])+'$dB$')
	# plot the polynomial representation of this exit curve 
	I_A_amp = np.linspace(0,1,101).reshape(-1,1)
	ones = np.ones((101, 1))
	# multiply the polynomial coefficients by the different values of I_A_amp and it's powers
	I_E_amp = np.sum(poly_coeff * np.concatenate((ones, I_A_amp, I_A_amp**2, I_A_amp**3),axis=1), axis=1)
	ax.plot(I_A_amp, I_E_amp, 'r--', label='Polynomial fit to EXIT chart')
	
	plt.xlabel('$I_A$')
	plt.ylabel('$I_E$')
	plt.legend(loc=6, prop={'size': 7})
	plt.title(str(poly_coeff))
	plt.show()
	#plt.savefig('polynomial_threshold_L256M32R1P4Bins350Threshold0_7_200reps_Plcorrect.pdf')



if __name__ == "__main__":
	t0=time.time()
	'''
	######
	# plotting the EXIT chart for the AMP decoder for a range of SNR
	sparcparams = SPARCParams(L=256, M=32, sigma=None, p=4, r=1, t=64)
	amp_exit_curve(sparcparams, low_snr_dB=10, high_snr_dB=13, repeats=1, threshold=0.7, x_axis_points=10, poly_curve=0, bin_number=350, import_data=False, export_csv_filename=None, import_csv_filename=None)
	'''
	
	###################################
	# Plotting the combined EXIT chart for amp and an LDPC with VND of different degrees
	# use this to see if the exit curves intersect
	R_ldpc = 1/2
	d_c=7
	# the average value of d_v
	dv_bar = (1-R_ldpc)*d_c
	
	# This is for M=64 L=256 R_sparc=1 snr=10dB
	poly_coeff = np.array([0.54368203, -0.3163521, 0.59025681, 0.23228106])
	# working with an ldpc code where d_v,1=2, d_v,2=4, d_v,3=18
	d_v1=2
	d_v2=7
	d_v3=12

	a_v = np.zeros(30)
	# have 2 planes which define the values of the a_v
	a, b = (d_v1, d_v2, d_v3, -dv_bar), (1, 1, 1, -1)
	# point on line of intersection point and direction of line d
	point, d = plane_intersect(a, b)

	a_v[d_v3] = 0.1
	A = np.array([[d_v1,d_v2],[1,1]])
	b = np.array([dv_bar-d_v3*a_v[d_v3], 1-a_v[d_v3]])
	# solve the simulataneous equations with a3 fixed
	print(A)
	print(b)
	x = np.linalg.solve(A,b)
	a_v[d_v1] = x[0]
	a_v[d_v2] = x[1]
	print(a_v)
	a_v = np.round(a_v, 3)
	
	#a_v[d_v1] = 29/40
	#a_v[d_v2] = 11/40
	#a_v[d_v3] = 0/40
	print(a_v)
	plot_amp_ldpc_exit(poly_coeff, a_v, d_c, R_ldpc)
	
	'''
	######## Looking for a feasible degree distribution
	a1, a2 = np.meshgrid(np.linspace(0,1,2),linspace(0,1,2))
	a3_1 = 1/d_v3 *(4/3 - d_v1*a1 - d_v2*a2)
	a3_2 = 1 - a1 - a2

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(a1,a2,a3_2)
	ax.plot_surface(a1,a2,a3_1)
	ax.set_xlabel("a1")
	ax.set_ylabel("a2")
	ax.set_zlabel("a3")
	plt.show()
	'''
	#plot_amp_ldpc_exit(poly_coeff, dv_count, a_v, d_c, R_ldpc)
	
	'''
	##### just plotting one set of histograms
	X = gen_bits(int(L*logm))
	print(X)
	EbN0_dB = 7
	snr = (10**(EbN0_dB/20))/(1/(2*R))
	snr_dB = 20*np.log10(snr)
	I_a = 0.8
	P=4
	snr = 10**(snr_dB/20)
	C = 0.5 * np.log2(1 + snr)
	sparcparams = SPARCParams(L=L, M=M, sigma=None, p=P, r=R_sparc, t=64)#, a=R_sparc/C, C=C, f=R_sparc/C)

	E = calc_E(X, I_a, snr_dB, sparcparams, threshold=0.6)
	np.clip(E, -55, 55, out=E)

	index_neg1 = np.where(X==-1)[0]
	print(np.ndarray.min(E))
	print(E[index_neg1])
	PE_pos, PE_neg, mean_pos, mean_neg, var_pos, var_neg, _ = hist_E(X, E, bin_number=350, max_bin=60, min_bin=-60, plot=True, snr_dB=snr_dB)
	print("positive mean: ", mean_pos)
	print("negative mean: ", mean_neg)
	print("positive variance: ", var_pos)
	print("negative variance: ", var_neg)
	
	

	'''
	'''
	#############################
	# plot histogram of the sectionwise probabilities for the correct position
	I_a = 0.5
	repeats = 20
	beta_plot = np.array([])
	for i in range(repeats):
		# Histograms of the section-wise probabilities produced by a particular I_a
		X = gen_bits(int(L*logm))
		# convert to LLRs and then to sectionwise probabilities
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
		beta_0 = bp2sp2(bitwise_a, L, M)

		# find the correct beta which corresponds to these X
		# Need to convert X, which contains -1 and +1s to bits
		X = (X-1)*-1/2
	    # convert the bits to indices for encoding
		X_indices = bits2indices(X, M)
		assert len(X_indices)==L
		# Find the correct beta corresponding to X
		β = np.zeros((L*M, 1))
		for l in range(L):
		    β[l*M + X_indices[l]] = 1#np.sqrt(n * Pl[l])
		# record the indices of the correct non-zero entries. 
		indices_correct = np.nonzero(β)[0]
		
		# get the sectionwise probabilities assigned to the correct non-zero entries for each of the L sections
		beta_plot = np.append(beta_plot, beta_0[indices_correct])

	print(shape(beta_plot))
	#print(beta_plot)
	bin_edges = np.linspace(0, 1, 100)
	plt.figure(1)
	plt.hist(beta_plot, bins=bin_edges, density = False)
	plt.title("Histogram of the section-wise probabilities for the correct non-zero entry in beta")
	plt.show()
	'''
	
	