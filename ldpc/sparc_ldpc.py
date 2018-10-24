import numpy as np
import math as math
from pylab import *
import py.ldpc as ldpc


# Code taken from Adam's python tutorial. 
# I have removed the power allocation code as I'm using a uniform power allocation
# Can get this from the notebook if I decide I need this. 


try:
    from pyfht import fht_inplace
except ImportError:
    import warnings
    warnings.warn("Using very slow Python version of fht, please install pyfht")
    def fht_inplace(x):
        N = len(x)
        i = N>>1
        while i:
            for k in range(0, N, 2*i):
                for j in range(k, i+k):
                    ij = i|j
                    temp = x[j]
                    x[j] += x[ij]
                    x[ij] = temp - x[ij]
            i = i >> 1


def sub_fht(n, m, seed=0, ordering=None):
    """
    Returns functions to compute the sub-sampled Walsh-Hadamard transform,
    i.e., operating with a wide rectangular matrix of random +/-1 entries.

    n: number of rows
    m: number of columns

    It is most efficient (but not required) for max(m+1,n+1) to be a power of 2.

    seed: determines choice of random matrix
    ordering: optional n-long array of row indices in [1, max(m,n)] to
              implement subsampling; generated by seed if not specified,
              but may be given to speed up subsequent runs on the same matrix.

    Returns (Ax, Ay, ordering):
        Ax(x): computes A.x (of length n), with x having length m
        Ay(y): computes A'.y (of length m), with y having length n
        ordering: the ordering in use, which may have been generated from seed
    """
    assert n > 0, "n must be positive"
    assert m > 0, "m must be positive"
    w = 2**int(np.ceil(np.log2(max(m+1, n+1))))

    if ordering is not None:
        assert ordering.shape == (n,)
    else:
        rng = np.random.RandomState(seed)
        idxs = np.arange(1, w, dtype=np.uint32)
        rng.shuffle(idxs)
        ordering = idxs[:n]

    def Ax(x):
        assert x.size == m, "x must be m long"
        y = np.zeros(w)
        y[w-m:] = x.reshape(m)
        fht_inplace(y)
        return y[ordering]

    def Ay(y):
        assert y.size == n, "input must be n long"
        x = np.zeros(w)
        x.flat[ordering] = y
        fht_inplace(x)
        return x[w-m:]

    return Ax, Ay, ordering

def block_sub_fht(n, m, l, seed=0, ordering=None):
    """
    As `sub_fht`, but computes in `l` blocks of size `n` by `m`, potentially
    offering substantial speed improvements.

    n: number of rows
    m: number of columns per block
    l: number of blocks

    It is most efficient (though not required) when max(m+1,n+1) is a power of 2.

    seed: determines choice of random matrix
    ordering: optional (l, n) shaped array of row indices in [1, max(m, n)] to
              implement subsampling; generated by seed if not specified, but
              may be given to speed up subsequent runs on the same matrix.

    Returns (Ax, Ay, ordering):
        Ax(x): computes A.x (of length n), with x having length l*m
        Ay(y): computes A'.y (of length l*m), with y having length n
        ordering: the ordering in use, which may have been generated from seed
    """
    assert n > 0, "n must be positive"
    assert m > 0, "m must be positive"
    assert l > 0, "l must be positive"

    if ordering is not None:
        assert ordering.shape == (l, n)
    else:
        w = 2**int(np.ceil(np.log2(max(m+1, n+1))))
        rng = np.random.RandomState(seed)
        ordering = np.empty((l, n), dtype=np.uint32)
        idxs = np.arange(1, w, dtype=np.uint32)
        for ll in range(l):
            rng.shuffle(idxs)
            ordering[ll] = idxs[:n]
        

    def Ax(x):
        assert x.size == l*m
        out = np.zeros(n)
        for ll in range(l):
            ax, ay, _ = sub_fht(n, m, ordering=ordering[ll])
            out += ax(x[ll*m:(ll+1)*m])
        return out

    def Ay(y):
        assert y.size == n
        out = np.empty(l*m)
        for ll in range(l):
            ax, ay, _ = sub_fht(n, m, ordering=ordering[ll])
            out[ll*m:(ll+1)*m] = ay(y)
        return out

    return Ax, Ay, ordering

#SPARC dictionary
#Returns 2 functions Ab and Az which compute Aβ and A^Tz respectively.
def sparc_transforms(L, M, n):
    Ax, Ay, _ = block_sub_fht(n, M, L, ordering=None)
    def Ab(b):
        return Ax(b).reshape(-1, 1) / np.sqrt(n)
    def Az(z):
        return Ay(z).reshape(-1, 1) / np.sqrt(n)
    return Ab, Az


# the amp algorithm
def amp(y, σ_n, Pl, L, M, T, Ab, Az):
    P = np.sum(Pl)
    n = y.size
    β = np.zeros((L*M, 1))
    z = y
    last_τ = 0
    
    for t in range(T):
        τ = np.sqrt(np.sum(z**2)/n)
        if τ == last_τ:
            return β
        last_τ = τ
        
        s = β + Az(z)
        rt_n_Pl = np.sqrt(n*Pl).repeat(M).reshape(-1, 1)
        u = s * rt_n_Pl / τ**2
        max_u = u.max()
        exps = np.exp(u - max_u)
        sums = exps.reshape(L, M).sum(axis=1).repeat(M).reshape(-1, 1)
        β = (rt_n_Pl * exps / sums).reshape(-1, 1)
        z = y - Ab(β) + (z/τ**2) * (P - np.sum(β**2) / n)
    
    return β

#########End of code from Python Tutorial###########

# SPARC parameters
# l is the number of sections in the overall SPARC
# m is the number of columns per section in the SPARC
# p is the channel signal power, sigma2 is the channel noise power
# r is the user data rate, r_pa the power allocation design rate
# t is the maximum number of AMP iterations permitted
class SPARCParams:
	def __init__(self, L, M, sigma, p, r, r_pa, t):
		self.L = L
		self.M = M
		self.sigma = sigma
		self.p = p
		self.r = r
		self.r_pa = r_pa
		self.t = t

# LDPC Params needed to set up one of Jossy's code
class LDPCParams:
    def __init__(self, standard, r_ldpc, z, ptype='A'):
        self.standard = standard
        self.r_ldpc = r_ldpc
        self.z = z
        self.ptype = ptype


# results of an LDPC simulation
class LDPCResult:
	def __init__(self, errors_after_amp_unprotected, errors_after_amp_protected, errors_after_amp_parity,
		errors_after_ldpc_protected, errors_after_ldpc_parity, errors_after_post_amp_unprotected, iters_amp,
		iters_ldpc, iters_post_amp, ldpc_success):
		self.errors_after_amp_unprotected = errors_after_amp_unprotected
		self.errors_after_amp_protected = errors_after_amp_protected
		self.errors_after_amp_parity = errors_after_amp_parity
		self.errors_after_ldpc_protected = errors_after_ldpc_protected
		self.errors_after_ldpc_parity = errors_after_ldpc_parity
		self.errors_after_post_amp_unprotected = errors_after_post_amp_unprotected
		self.iters_amp = iters_amp
		self.iters_ldpc = iters_ldpc
		self.iters_post_amp = iters_post_amp
		self.ldpc_success = ldpc_success



#function to convert sectionwise probabilities to bitwise probabilities
# note: my sp2bp takes the entire beta vector, whereas adams code works only on 1 section. Change this if neccessary
def sp2bp(β, L, M: "must be power of 2"):
    #take in the section posterior probabilities β.
    #and the no. sections L and the size of each section M. 
    
    #initialise numpy array of zeros for all of the bitwise posteriors
    p = np.zeros(int(np.log2(M)*L))
    #loop through the L sections.
    for a in range(L):
        β_l = β[a*M:((a+1)*M)]
        # normalizing before passing in, so below not needed
        #c = np.sum(β_l) #calculate normalization constant for each section
        for logi in range(int(np.log2(M))):
            b = int((a+1)*np.log2(M) - logi - 1)
            i = 2**logi
            k = i
            while k<M:
                #note shift in range of j due to βl being an array which is indexed from 0 not 1. 
                for j in range(k, k+i):
                    p[b] = p[b] + β_l[j]
                k = k + 2*i
    return p
    


# assumes that the left most bit (lmb) is the most significant bit (msb)
# note that an input of all zero bits for one section will give an output of 0
# So the first column is being indexed by 0
def bits2indices(bits, m: "m must be a power of 2"):
	
    logm = int(math.log(m,2))
    # assert that the bits can be split exactly into sections of size logm
    assert len(bits)%logm==0
    # find number of sections required for the bits
    sections = int(len(bits)/logm)
    # set up array to store the indices
    indices = []

    for i in range(sections):
        # get the bits corresponding to the current section
        section_bits = bits[i*logm:(i+1)*logm]
        # variable to store the value of the bits in indice
        digit = 0
        for j in range(len(section_bits)):
            if section_bits[j]:
                digit += 2**(len(section_bits)-j-1) 
        indices.append(digit)

    return indices

def amp_ldpc_sim(sparcparams: SPARCParams, ldpcparams: LDPCParams = None):
    #Get the SPARC parameters from the struct
    L = sparcparams.L
    M = sparcparams.M
    P = sparcparams.p
    sigma = sparcparams.sigma
    r_sparc = sparcparams.r
    r_pa_sparc = sparcparams.r_pa
    T = sparcparams.t

    # calculate some additional parameters 
    # Compute the SNR, capacity, and n, from the input parameters
    snr = P / sigma**2
    #C = 0.5 * np.log2(1 + snr)
    n = int(L*np.log2(M) / r_sparc)
    # compute additional parameters
    logm = np.log2(M)
    total_bits = int(logm*L)
    
    # Generate the power allocation
    # Pl = pa_iterative(L, L, sigma, P, R_PA)
    # uniform power allocation across sections
    Pl = P/L * np.ones(L)


    if ldpcparams == None:
        # set params so that no ldpc encoding or decoding is used
        nl = 0
        kl = 0
        ldpc_bits = []
    else: 
        standard = ldpcparams.standard
        r_ldpc = ldpcparams.r_ldpc
        z = ldpcparams.z
        ptype = ldpcparams.ptype
        # initialise the ldpc code
        ldpc_code = ldpc.code(standard, r_ldpc, z, ptype)
        nl = ldpc_code.N
        kl = ldpc_code.K

        assert nl<=(L*logm)
        # we want the ldpc to cover a complete number of sections and not just part of the final section it covers. 
        assert nl%logm == 0
        #(L, M, sigma, P, R, T, R_PA, R_LDPC):

        # Generate random message in 2 stages
        # First generate ldpc bits
        protected_bits = np.random.randint(0, 2, kl).tolist()
        assert len(protected_bits) == kl

        # Encode the protected_bits using the ldpc code 
        ldpc_bits = ldpc_code.encode(protected_bits).tolist() 

    # Generate the remaining bits required
    unprotected_bits = np.random.randint(0, 2, int(L*logm-nl)).tolist()

    # concatenate the ldpc and unprotected bits      
    sparc_bits = unprotected_bits+ldpc_bits

    assert len(sparc_bits)==total_bits

    # convert the bits to indices for encoding
    sparc_indices = bits2indices(sparc_bits, M)

    assert len(sparc_indices)==L
       
    # Generate the SPARC transform functions A.beta and A'.z
    Ab, Az = sparc_transforms(L, M, n)
    
    # Generate our transmitted signal X
    β_0 = np.zeros((L*M, 1))
    for l in range(L):
        β_0[l*M + sparc_indices[l]] = np.sqrt(n * Pl[l])
    x = Ab(β_0)
    
    # check that the power has been allocated uniformly. This should be approx equal to the snr. 
    print("Average power/snr is ", np.mean(x**2)/snr)

    # Generate random channel noise and then received signal y
    z = np.random.randn(n, 1) * sigma
    y = (x + z).reshape(-1, 1)
        
    # Run AMP decoding
    β = amp(y, sigma, Pl, L, M, T, Ab, Az).reshape(-1)
    # Output the results of just AMP decoding
    
    # Keep below so I can compare it to using the outer ldpc code. Actually below doesn;t work anymore as I generated bits initially. 
    # Convert decoded beta back to a message
    rx_message = []
    for l in range(L):
        idx = np.argmax(β[l*M:(l+1)*M])
        rx_message.append(idx)
    
    # Compute fraction of sections decoded correctly with just AMP decoding
    correct_amp = np.sum(np.array(rx_message) == np.array(sparc_indices)) / L
    print("Fraction of sections decoded correctly with amp: ", correct_amp)
    
    # Compute BER note: np.sum will be deprecated soon, so look up alternative if this code stops working. 
    # DeprecationWarning: Calling np.sum(generator) is deprecated, and in the future will give a different result. Use np.sum(np.from_iter(generator)) or the python sum builtin instead.
    # If this function starts acting weird, look into the above. Or when I transfer this function to my proper code, try to fix this. 
    ber_amp = sum(bin(a^b).count('1') for (a, b) in zip(sparc_indices, rx_message))/(L*logm)


    if ldpcparams==None:
        ber_ldpc = None
    else:
        # Get the sectionwise posterior probabilities by dividing β by the power in each section. 
        # This converts each section in β to a valid probability distribution.
        sectionwise_posterior = β/np.sqrt(n*P/L)


        # convert sectionwise to bitwise posterior probabilities
        bitwise_posterior = sp2bp(sectionwise_posterior, L, M) 
        #print(np.sum(bitwise_posterior))


        # only want to use the bitwise posterior for the sections we've applied the ldpc code to
        # recenter the bitwise posterior around zero for the decoding
        LLR = np.log(1-bitwise_posterior[int(total_bits-nl):])- np.log(bitwise_posterior[int(total_bits-nl):])

        (app, it) = ldpc_code.decode(LLR)
        #print((np.ones(nl)/2)-bitwise_posterior[int((L*logm)-nl):])
        test_output = (bitwise_posterior[int(total_bits-nl):]>0.5)
        print("Incorrect bits before ldpc decoding: ", np.sum(test_output!=ldpc_bits))

        # check that the inequality is the right way round in below
        v_output = (app<0.0)

        print("Incorrect bits in ldpc decoding: ", np.sum(v_output!=ldpc_bits))

        # convert the bits to indices
        beta_output_ldpc = bits2indices(v_output, M)

        beta_output = rx_message
        beta_output[int(L - nl/logm):] = beta_output_ldpc
        print(np.sum(rx_message[int(L - nl/logm):]!= beta_output_ldpc))
        # compute fraction of sections decoded correctly with AMP followed by SPARC decoding
        correct_ldpc = np.sum(np.array(beta_output) == np.array(sparc_indices))/L
        print("Fractions of sections decoded correctly with amp and ldpc: ", correct_ldpc)

            
        # compute BER for amp and ldpc
        ber_ldpc = sum(bin(a^b).count('1') for (a, b) in zip(sparc_indices, beta_output))/(L*np.log2(M))

    # overall rate of the code
    R = (L*logm - (nl-kl))/n

    #Compute Eb/N0
    EbN0 = 1/(2*R) * (P/sigma**2)
   
    return ber_amp, ber_ldpc, R




if __name__ == "__main__":
    '''c = np.array([0.2, 0.3, 0.1, 0.4, 0.6, 0.4, 0, 0])
    print(sp2bp(c, 2, 4))

    '''
    M = 512
    L = 512
    sigma = 1
    P = 4
    r_sparc = 1
    r_pa_sparc = 1
    T = 64
    sparcparams = SPARCParams(L, M, sigma, P, r_sparc, r_pa_sparc, T)
    
    
    # The ldpc sections must be less than L and must be divisible by 8 to give an integer value for z
    ldpc_sections = 400
    # nl is the number of bits that will be LDPC bits. 
    # This is given by the no. bits per sections * the no. sections
    # Note the no. sections must be at least 8 to give something divisible by 24
    nl = np.log2(512) * ldpc_sections
    # Can get from nl to z. z determines the size of the protograph matrix
    z = int(nl/24)
    standard = '802.16'
    rate = '1/2'
    ldpcparams = LDPCParams(standard, rate, z)

    (ber_amp, ber_ldpc, R) = amp_ldpc_sim(sparcparams, ldpcparams)
    print("ber_amp: ", ber_amp)
    print("ber_ldpc: ", ber_ldpc)
    print("R: ", R)


    '''
    #Code to plot the performance of the sparc code on it's own. 
    # fix rate
    R = 1

    σ_n = 1.0
    # shannon limit on snr. i.e. snr at capacity
    snrc = 2**(2*R)-1
    
    # Unsure what values to choose for T and R_PA. Just copied previous example.
    T = 64
    R_PA = R
    # store the error rate for each snr
    error_rate = []
    repeats = 1
    for snr in range(snrc,10):
        # calculate the new value of P to adjust the snr
        P = snr * σ_n**2
        error = []
        sparcparams = SPARCParams(512, 512, 1, P, 1, 1, 64)
        for i in range(repeats):
            ber = amp_ldpc_sim(sparcparams)[0]
            error.append(ber)
        aver_error = np.sum(error)
        error_rate.append(aver_error/repeats)
        
    plot(range(snrc,10), error_rate)
    xlabel('snr')
    title("Parameters: L=512, M=512, σ_n=1, R=1, T=64, R_PA=1")
    ylabel('BER')
    show()'''


