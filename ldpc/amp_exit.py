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
	def __init__(self, X, E, I_a, SNR_dB):
		self.X = X
		self.E = E
		self.I_a = I_a 
		self.SNR_dB = SNR_dB
	def __eq__(self, other): 
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
			# set the position corresponding to the probability greater than the threshold to 1
			beta_0[(l*M)+idx[0]] = 1	
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

# Function produces y from X a series of bits for inputting into the amp() function
# X: array of LlogM input bits
# L, M, N are all SPARC parameters
# sigma_w is the variance of the noise
# P is the power that will be uniformly allocated
# returns: y, Ab, Az, and Pl
def prep_y(X, L, M, n, sigma_w, P, a=None, f=None, C=None):
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
def calc_E(X, I_a, snr_dB, sparcparams, csv_filename=None, threshold=0.5):
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
	# Generate y, Ab, and Az to pass into amp()
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
		# clip the bitwise_e to avoid a divide by zero error in the log
		#np.clip(bitwise_e, 0.00000000000000000000000001, 1-0.00000000000000000000000001, out=bitwise_e)
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

# given E and X return the histograms of E
# X = generated input bits
# bin_number, max_bin, min_bin determine the number of histogram bins and the range of the bins
# plot is a bool. Set true to plot histograms
# Return: histograms of E, and their means and variances. 
def hist_E(X, E, bin_number=500, max_bin=40, min_bin=-40, plot=False, snr_dB='Not given'):
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
		#plt.title("Normalised histogram of $P(E(X_i)|X_i=+1)$."+" $SNR=$"+str(snr_dB)+"dB")

		plt.figure(2)
		plt.hist(E[index_neg1], bins=bin_edges, density = True)
		#plt.title("Normalised histogram of $P(E(X_i)|X_i=-1)$."+" $SNR=$"+str(snr_dB)+"dB")
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

# Parameters: I_a array containing all the I_a for the EXIT curve. 
# I_e array containing the corresponding I_e values. 
# Returns the 3rd order polynomial that best fits an EXIT curve based on LS
# i.e. returns array of the polynomial coefficients. 
# return has order c0,c1,c2,c3 where c0 is the constant and c3 is the coefficient of I_a^3
def polynomial(I_a, I_e):
	a = np.zeros((len(I_a), 4))
	for i in range(4):
		a[:,i] = I_a**i
	c = np.linalg.lstsq(a, I_e, rcond=-1)[0]
	return c	

# returns I_E_VND, the output from the combined AMP decoder and VND
# The amp exit curve is represented by a 3rd order polynomial equation
# with the coefficients passed in as array poly_coeff
# d_v is the degree of the variable node in question
def I_E_VND_amp(I_A_VND, d_v, poly_coeff):
	I_A_amp = J(np.sqrt(d_v)*J_inverse(I_A_VND))
	I_E_amp = np.sum(poly_coeff * np.array([1, I_A_amp, I_A_amp**2, I_A_amp**3]))
	# clip I_E_amp at a max of I_E_amp as it cannot exceed this
	I_E_amp = np.clip(I_E_amp, a_min=None, a_max=0.9999)
	I_E_VND = J(np.sqrt((d_v-1)*(J_inverse(I_A_VND))**2+(J_inverse(I_E_amp))**2))
	return I_E_VND

# calculate the array I_E_VND as a function of I_A_VND and I_A_amp over a range of d_vs.
# a_v the fraction of variable nodes with degree equal to its index.
# b_v the fraction of edges incident to variables nodes of degree equivalent to the index
# poly_coeff the 3rd order polynomial fit to the amp exit chart.
def I_E_VND_amp_array(I_A_VND_array, a_v, b_v, poly_coeff):
	I_E_VND = np.zeros(I_A_VND_array.shape)
	for i in range(len(a_v)):
		# a_v gives the fraction of variable nodes with degree equal to its index.
		# when this is zero means there are no v nodes with this degree so skip it
		if a_v[i]!=0:
			for j in range(len(I_E_VND)):
				I_E_VND[j] = I_E_VND[j]+b_v[i]*I_E_VND_amp(I_A_VND_array[j], i, poly_coeff)
	return I_E_VND

# Plot the exit charts for the combined amp and VND and the CND
# This function is for a fixed value of CND.
# the ldpc curve can have variable nodes of different degrees
# Params:
# poly_coeff - the coefficients of the 3rd order polynomial representating the amp exit curve
# a_v array giving the fraction of variable nodes of the degree equal to its index 
# R_ldpc the rate of the ldpc code
def plot_amp_ldpc_exit(poly_coeff, a_v, d_c, R_ldpc):
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
	plt.plot(I_A_VND, I_E_VND, color='k', linestyle='-', label = 'VND')
	plt.plot(I_E_CND, I_A_CND, color='r', linestyle='-', label = 'CND')
	plt.xlabel('$I_{A,VND}, I_{E,CND}$', fontsize=18) 
	plt.ylabel('$I_{E,VND}, I_{A,CND}$', fontsize=18)
	plt.tight_layout()
	plt.legend()
	plt.show()

# plot the exit charts for the combined amp and VND and the CND 
# Both the VND and CND can have variable degree distributions. 
#def plot_amp_ldpc_exit_varCND(poly_coeff, a_v, )

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


if __name__ == "__main__":
	t0=time.time()
	L=256
	M=32
	logm = np.log2(M)
	R_sparc = 1
	export = True

	'''
	###################################
	# Plotting the combined EXIT chart for amp and an LDPC with VND of different degrees
	R_ldpc = 3/8
	d_c=6
	# the average value of d_v
	dv_bar = (1-R_ldpc)*d_c
	
	# This is for M=64 L=256 R_sparc=1 snr=10dB
	poly_coeff = np.array([0.43764836, 0.71227327, -2.57287966, 2.4402426])
	# working with an ldpc code where d_v,1=2, d_v,2=4, d_v,3=18
	d_v1=2
	d_v2=5
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
	'''
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
	#plot_amp_ldpc_exit(poly_coeff, dv_count, a_v, d_c, R_ldpc)
	'''	
	'''
	# just plotting one set of histograms
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
	
	# plotting the EXIT chart for the AMP decoder for a range of SNR
	bin_number = 350
	repeats = 200
	datapoints = 3
	x_axis_points=30
	I_a_range = np.linspace(0, 0.99, x_axis_points)
	P = 4
	if export==False:	# if not exporting, want to import E into a dictionary
		# sometimes need to set the repeats in here higher than the actual repeats to ensure all the data is imported. I don't understand why!
		imported_E_dict	= import_E_fromfile(fileName = 'exit_charts/E_data__L512_M512_20reps_500bins_r_1.csv', datapoints = datapoints, repeats = repeats, Llogm = int(L*logm))
		print(len(imported_E_dict))
	# accumulative values of I_e for each snr value
	I_e_accum = np.zeros((datapoints,x_axis_points))
	#EbN0_dB = np.linspace(5, 8, datapoints)
	#snr = (10**(EbN0_dB/20))/(1/(2*R))
	# work in snr for the EXIT charts as then don't have to work about EbN0 and rate. 
	snr_dB = np.linspace(10, 12, datapoints)
	for k in range(repeats):
		j=0
		for s_dB in snr_dB:
			# channel capacity
			snr = 10**(s_dB/20)
			C = 0.5 * np.log2(1 + snr)
			sparcparams = SPARCParams(L=L, M=M, sigma=None, p=P, r=R_sparc, t=64)#, a=R_sparc/C, C=C, f=R_sparc/C)
			I_e = np.zeros(x_axis_points)
			i=0
			for I_a in I_a_range:
				a = None
				if export==True: # if export if false, we already have these values.
					# randomly generate bits
					X = gen_bits(int(L*logm))
					#print(X)

					# generate the histograms for E and some statistics about them
					E = calc_E(X, I_a, s_dB, sparcparams, threshold=0.85)#, csv_filename='E_data_hardinit_LM512R1P4Bins125Threshold0_45.csv')
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
				PE_pos, PE_neg, mean_pos, mean_neg, var_pos, var_neg, bin_width = hist_E(X, E, bin_number=bin_number, max_bin=60, min_bin=-60)#, plot=True, snr_dB=s_dB)
	
				# calculate I_e
				I_e[i] = calc_I_e(PE_pos, PE_neg, bin_width)
				i=i+1
			I_e_accum[j,:] = I_e_accum[j,:] + I_e
			j=j+1
	I_e_accum = I_e_accum/repeats
	# calculate the polynomial coefficients for the 3rd order polynomial representation of the exit curve
	curve=1
	poly_coeff = polynomial(I_a_range, I_e_accum[curve,:]) 
	# Note that this will give the constant first and the coefficient of I_a**3 last
	print("The coefficients for the polynomial are: ",poly_coeff)
	print("Wall clock time elapsed: ", time.time()-t0)

	plt.figure(1)
	fig, ax = plt.subplots()
	ax.plot(I_a_range, I_e_accum[0,:], 'b--', label='$SNR$='+str(snr_dB[0])+'$dB$')	
	ax.plot(I_a_range, I_e_accum[1,:], 'k--', label='$SNR$='+str(snr_dB[1])+'$dB$')
	ax.plot(I_a_range, I_e_accum[2,:], 'm--', label='$SNR$='+str(snr_dB[2])+'$dB$')
	#ax.plot(I_a_range, I_e_accum[3,:], 'c--', label='$SNR$='+str(snr_dB[3])+'$dB$')
	#ax.plot(I_a_range, I_e_accum[4,:], 'r--', label='$SNR$='+str(snr_dB[4])+'$dB$')
	plt.xlabel('$I_A$')
	plt.ylabel('$I_E$')
	plt.legend(loc=6, prop={'size': 7})
	#plt.title("The EXIT chart for the AMP decoder")
	#plt.savefig('amp_exitchart_L128_M4_40reps_500bins_r1_5_P2.png')	
	plt.savefig('amp_exit_threshold_L256M32R1P4Bins350Threshold0_85_200reps.pdf')	
	#plt.show()
	
	plt.figure(2)
	fig, ax = plt.subplots()
	ax.plot(I_a_range, I_e_accum[curve,:], 'b--', label='$SNR$='+str(snr_dB[curve])+'$dB$')
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
	#plt.show()
	plt.savefig('polynomial_threshold_L256M32R1P4Bins350Threshold0_85_200reps.pdf')
	
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
	
	