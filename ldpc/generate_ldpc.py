import numpy as np

class ParityMatrix:
	def __init__(self, dv_degrees, a_v, d_c):
		# dv_degrees contains the variable node degree distribution. e.g. dv_degrees=[2, 7, 12]
		self.dv_degrees = dv_degrees
		# contains the fraction of each degree in dv_degrees e.g a_v=[0.8, 0.1, 0.1]
		self.a_v = a_v
		# the degree of the check node. The code is check regular so this can only take 1 value
		self.d_c = d_c
		# the average variable node degree
		self.dv_bar = self.dv_bar()
		return

	def dv_bar(self):
		dv_bar=0
		for i in range(len(self.a_v)):
			dv_bar = dv_bar + self.a_v[i]*self.dv_degrees[i]
		# rounding to 4dp due to small errors that can be introduced in equation above
		# later causing K to not be an integer
		return np.round(dv_bar, decimals=4)

	# N and K are the dimensions of the LDPC code
	# This function calculates the parity check matrix 
	def generate(self, N, K):
		# the CND distribution is check regular. Use this to keep track of how many 1s you
		# still need in each row 
		dcnd = np.ones(N-K)*self.d_c
		# randomly sample the VND distribution from a_v
		dvnd = np.random.choice(self.dv_degrees, size=N, p=a_v)
		print("dvnd is: ", dvnd)
		# if there aren't enough rows to satisfy some of the required degrees
		if N-K<np.amax(self.dv_degrees):
			print("Please choose an N-K greater than ", np.amax(self.dv_degrees))
			return

		parity_matrix = np.zeros((N-K, N))
		# loop through the columns of the parity check matrix
		for col in range(N):
			#print(col)
			dv = dvnd[col]
			#print(dcnd)
			if np.sum(dcnd)<dv or np.count_nonzero(dcnd)<dv:
				# if there aren't enough check node edges left to create all the required columns
				# or if there the number of nonzero entries if less that dv (i.e. which would mean
				# we couldn't sample dv things with the replace=False) 
				# add one to every entry of dcnd. This means we can continue to add columns 
				# but they may have more than d_c entries. By only adding one when required, 
				# we hopefully avoid getting many more than d_c entries in any particular row
				dcnd = dcnd + 1
				#print(dcnd)
			# This should pick random rows of the parity check matrix to be nonzero
			# by setting replace to false, you should never pick the same row twice. 
			idx_ones = np.random.choice(np.arange(0,N-K),size=dv,replace=False,p=dcnd/np.sum(dcnd))

			# decrement all the check node distribution that have had a variable node connected. 
			dcnd[idx_ones] += -1 
			parity_matrix[idx_ones,col] = 1

		return parity_matrix.astype(int)

		# returns a protograph matrix which when used with z (also returned) will give 
		# a parity check matrix of dimensions N and K
		# N and K must be divisible by 3, 4 or 5 and must desirably not be too small to ensure
		# all dv degrees get included in the parity check matrix. 
	def generate_protograph(self, N, K):
		# to generate a protograph we must divide N and N-K by 3 or more so that
		# producing the desired parity check matrix involves using a z of 3 or greater
		# where z is the dimension of the shifted identity or zero matrix that replaces 
		# each entry in a parity check matrix. 
		# n and n-k are the dimensions of the protograph matrix 
		# z is the required dimension of the shifted/zero matrices to give Nx(N-K) parity matrix
		if N%3==0 and (N-K)%3==0:
			n=N/3
			k=K/3
			z=3
		elif N%4==0 and (N-K)%4==0:
			n=N/4
			k=K/4
			z=4
		elif N%5==0 and (N-K)%5==0:
			n=N/5
			k=K/5
			z=5
		else:
			print("choose a different value of N and K")
			return 
		parity = self.generate(int(n), int(k))
		# replace all zeros by -1, i.e. the all zero matrix. 
		np.place(parity, parity==0, -1)
		#print(parity)
		# find the location of all the ones in the parity matrix
		idx_ones = np.where(parity>0)
		# count the number of ones. Note idx_ones contains 2 arrays of length total_ones
		total_ones = len(idx_ones[0])
		# randomly generate the protograph permatations. I have decided to max this at 100.
		proto_shift = np.random.choice(np.arange(0,100), size=total_ones)
		# replace all ones by a random permatation to produce the protograph matrix
		parity[idx_ones]=proto_shift
		return parity.astype(int), z


##################### Testing generate() and generate_protograph()
dv_degrees = np.array([2, 7, 12])
a_v = np.array([0.8, 0.1, 0.1])
d_c = 7
pm = ParityMatrix(dv_degrees, a_v, d_c)
dv_bar=pm.dv_bar
N=1280/8
print(dv_bar)
K=N-(dv_bar*N/d_c)
print(K)
#print(pm.generate(N, int(K)))
np.set_printoptions(threshold=np.nan)
print(pm.generate_protograph(N,int(K)))

