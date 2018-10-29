# This contains a series of functions that have been removed from other files
# but I don't want to delete in case they become useful in the future

############## Functions from sparc_ldpc.py ################
# performs sparc encoding and decoding without using an outer ldpc
# deprecated as my ldpc_amp_sim should do everything this does when ldpc_params=None
def amp_sim(L, M, σ_n, P, R, T, R_PA):
    # Compute the SNR, capacity, and n, from the input parameters
    #snr = P / σ_n**2 don't actually need snr and C
    #C = 0.5 * np.log2(1 + snr)
    # number of channel uses by the sparc code. 
    n = int(L*np.log2(M) / R)
    
    # Generate the power allocation
    # Pl = pa_iterative(L, L, σ_n, P, R_PA)
    # uniform power allocation across sections
    Pl = P/L * np.ones(L)
    
    # Generate random message in [0..M)^L
    tx_message = np.random.randint(0, M, L).tolist()
    
    # Generate the SPARC transform functions A.beta and A'.z
    Ab, Az = sparc_transforms(L, M, n)
    
    # Generate our transmitted signal X
    β_0 = np.zeros((L*M, 1))
    for l in range(L):
        β_0[l*M + tx_message[l]] = np.sqrt(n * Pl[l])
    x = Ab(β_0)

    # check that the power has been allocated uniformly. 
    print(np.mean(x**2))
    
    # Generate random channel noise and then received signal y
    z = np.random.randn(n, 1) * σ_n
    y = (x + z).reshape(-1, 1)
        
    # Run AMP decoding
    β = amp(y, σ_n, Pl, L, M, T, Ab, Az).reshape(-1)
    #The above is the section wise probabilities. I need to convert these to bitwise probabilities. 
    
    
    #print(bitwise_posterior(β, L, M))
    
    # Convert decoded beta back to a message
    rx_message = []
    for l in range(L):
        idx = np.argmax(β[l*M:(l+1)*M])
        rx_message.append(idx)
    
    # Compute fraction of sections decoded correctly
    correct = np.sum(np.array(rx_message) == np.array(tx_message)) / L
    
    # Compute BER note: np.sum will be deprecated soon, so look up alternative if this code stops working. 
    ber = np.sum(bin(a^b).count('1')
                 for (a, b) in zip(tx_message, rx_message))/(L*np.log2(M))
    
    # Compute Eb/N0
    EbN0 = 1/(2*R) * (P/σ_n**2)
    
    return ber
    '''return {
        "L": L, "M": M, "sigma_n": σ_n, "P": P, "R": R, "T": T,
        "snr": snr, "C": C, "n": n, "fc": correct, "R_PA": R_PA,
        "ber": ber, "EbN0": EbN0, "ser": 1-correct
    }'''


# copying over from Adams rust. In the end it made more sense to adapt the python code he'd already written 
#Fn to run point-to-point SPARC-AMP simulation with an LDPC outer code protecting the flat sections.
def lpdc_sim(sparcparams: SPARCParams, seed)	-> LDPCResult:	
	#Get the SPARC parameters from the struct
	L = sparcparms.L
	M = sparcparams.M
	p = sparcparams.p
	sigma2 = sparcparams.sigma2
	r = sparcparams.r
	r_pa = sparcparams.r_pa
	t = sparcparams.t


	# pick ldpc code, get it's parameters


	#compute remaining parameters
	logm = int(log(m,2))
	n = int((l * logm - (nl-kl))/r) #nl and kl are the ldpc parameters n and k


	#set up random numbers generator
	np.random.seed(seed)

	#generate random user bits to transmit
	num_sparc_bits = l * logm
	num_user_bits = num_sparc_bits - (nl-kl)
	user_bits = np.random.randint(0, 2, num_user_bits).tolist()

	#split bits into amp bits and those that will be LDPC encoded
	#What are the AMP bits? Aren't all bits AMP decoded.
	amp_bits = user_bits[:(num_user_bits-kl)]
	ldpc_user_bits = user_bits[(num_user_bits-kl):]

	#LPDC encode those bits, via a mapping into bytes and back (this will determine if I need to bytes2bits fns etc.)

	#gather the AMP bits and the LDPC codeword bits (the "SPARC" bits)
	sparc_bits = amp_bits + ldpc_code_bits
	#turn all the SPARC bits into column choices
	#This is when I will need bits2indices
	sparc_indices = bits2indices(sparc_bits, m)


	#Make PA and FHT
	#PA is the power allocation? And FHT is the design matrix?
	#What are these? They are used for next comment. Is one the design matrix?

	#Encode sparc indecs into an actual codeword. 

	#generate channel noise and sum. 

	#Run AMP decoded

	#Gather the hard decision indices

	#Convert indices back to bits

	#Count bit errors in AMP alone. 

	#Convert all column probabilities into bit probabilities.

	#grab the LDPC bit probs, convert to LLRs, decode LDPC codeword

	#Count bit error in LDPC section

	#If LDPC worked, and there are errors in the unprotected bits, try running AMP again

		#Make a new beta_ldpc that only has the ldpc sections set

		#Compute a y that already has all the successfully decoded ldpc sections removed

		#Run AMP over the non-LDPC sections from scratch

		#Gather hard decision indices

	#Convert indices back to bits

	#Count bit errors in AMP unprotected bits

	#RETURN the LDPC results using the above defined struct
