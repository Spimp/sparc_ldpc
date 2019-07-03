import numpy as np
import math as math
from pylab import * 
import py.ldpc as ldpc
import matplotlib.pyplot as plt
import time
import csv


# Need to import the protograph matrix from ldpc. Or copy it in
# 0 corresponds to the identity matrix and -1 to the all zero matrix
proto = np.array([
                        [17, 13, 8, 21, 9, 3, 18, 12, 10, 0, 4, 15, 19, 2, 5, 10, 26, 19, 13, 13, 1, 0, -1, -1],
                        [3, 12, 11, 14, 11, 25, 5, 18, 0, 9, 2, 26, 26, 10, 24, 7, 14, 20, 4, 2, -1, 0, 0, -1],
                        [22, 16, 4, 3, 10, 21, 12, 5, 21, 14, 19, 5, -1, 8, 5, 18, 11, 5, 5, 15, 0, -1, 0, 0],
                        [7, 7, 14, 14, 4, 16, 16, 24, 24, 10, 1, 7, 15, 6, 10, 26, 8, 18, 21, 14, 1, -1, -1, 0]
                    ])

# J and J_inverse based on approximation in the appendix of Design of LDPC Codes for modulation and detection.
def J_inverse(I):
	assert I>=0 and I<=1
	if (I<=0.3646):
		return 1.09542*(I**2) + 0.214217*(I) + 2.33727*np.sqrt(I)
	else:
		return -0.706692*np.log(0.386013*(1-I)) + 1.75017*(I) 

def J(sigma):
	assert sigma>=0
	if (sigma<=1.6363):
		return -0.0421061*(sigma**3) + 0.209252*(sigma**2) + -0.00640081*(sigma)
	elif (sigma<10):
		return 1- np.exp(0.00181491*(sigma**3) - 0.142675*(sigma**2) - 0.0822054*sigma + 0.0549608)
	else:
		return 1

def I_E_VND(I_A, dv, EbN0, R):
	sigma2_ch = 8*R*EbN0
	sigma = np.sqrt((dv-1)*((J_inverse(I_A))**2) + sigma2_ch)
	return J(sigma)

def I_A_CND(I_E, dc):
	sigma = (J_inverse(1-I_E))/np.sqrt(dc-1)
	return 1-J(sigma)

# I_E_CND=1-I_E_REP So this function is used to calculate I_E_CND
def I_E_REP(I_A, dc):
	sigma = np.sqrt(dc-1)*J_inverse(1-I_A)
	return J(sigma)


if __name__ == "__main__":
	# as often dealing with fractions, the value of z shouldn't matter
	z=1
	# rate of the ldpc code
	R=5/6
	EbN0_dB = 7
	EbN0 = 10**(EbN0_dB/20)

	num_rows = proto.shape[0]
	num_cols = proto.shape[1]
	# Find the degree of the variable nodes and store in array
	# There are 5 different degrees you can get
	# counts the number of variable nodes of each degree
	dv_count = np.zeros(num_rows+1)
	# gives the degree of the variable nodes
	dv_degrees = linspace(0, 4, 5)
	for col in proto.T:
		# add 1 to col so we count the number of '-1' in proto. 
		dv_count[np.count_nonzero(col+1)] += 1

	# Find degree of check nodes and store in array
	# There are 25 different degrees you could get so
	# counts the number of check nodes of each degree
	dc_count = np.zeros(num_cols+1)
	# gives the degree of the check nodes
	dc_degrees = linspace(0,24,25)
	for row in proto:
		# add 1 to row so we count the number of '-1' in proto.
		dc_count[np.count_nonzero(row+1)] += 1
	# Find average check node degree
	dc_aver = (1/num_rows)*sum(dc_count*dc_degrees)

	# Find b_i for the check nodes
	b_c = dc_count*dc_degrees/(num_rows*dc_aver)
	print(sum(b_c))
	# Find b_i for the variable nodes
	a_v = dv_count/sum(dv_count)
	dv_aver = sum(a_v*dv_degrees)
	b_v = (dv_degrees*a_v)/dv_aver  #((1-R)*dc_aver)
	print(sum(b_v))
	print(sum(a_v))

	# plot I_E,VND as a func of I_A
	I_Av = linspace(0, 1, 21)
	I_Ev = np.zeros(I_Av.shape)
	for i in range(len(dv_count)):
		# dv_count gives the num of variable nodes with degree equal to its index.
		# when this is zero means there are no v nodes with this degree so skip it
		if dv_count[i]!=0:
			for j in range(len(I_Ev)):
				I_Ev[j] = I_Ev[j]+b_v[i]*I_E_VND(I_Av[j], i, EbN0, R)


	# plot I_A,CND as a func of I_E and plot this on the same graph
	'''I_Ec = linspace(0.1, 1, 9)
	I_Ac = np.zeros(I_Ec.shape)
	for k in range(len(dc_count)):
		if dc_count[k]!=0:
			for j in range(len(I_Ac)):
				I_Ac[j] = I_Ac[j] + b_c[k]*I_A_CND(I_Ec[j], k)'''

	# Find I_E,CND as a func of I_A and then plot this on the same graph 
	I_Ac = linspace(0.1, 1, 19)
	I_Ec = np.zeros(I_Ac.shape)
	for k in range(len(dc_count)):
		if dc_count[k]!=0:
			for j in range(len(I_Ec)):
				I_Ec[j] = I_Ec[j] + b_c[k]*I_E_REP(I_Ac[j], k)
	# subtract from 1 to complete the approximation. See equation 4 in Jossy's paper.
	I_Ec = 1-I_Ec

	fig, ax = plt.subplots()
	plt.plot(I_Av, I_Ev, color='k', linestyle='-', label = 'VND')
	plt.plot(I_Ec, I_Ac, color='r', linestyle='-', label = 'CND')
	plt.xlabel('$I_{A,VND}, I_{E,CND}$', fontsize=18) # at some point need to work out how to write this so it outputs properly
	plt.ylabel('$I_{E,VND}, I_{A,CND}$', fontsize=18)
	plt.tight_layout()
	plt.legend()
	plt.show()
    #plt.savefig('EXITchart_test.png')
	











