import numpy as np
import math as math
from pylab import * 
import py.ldpc as ldpc
import matplotlib.pyplot as plt
import time
import csv
from bitarray import bitarray


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
    Ax, Ay, ordering = block_sub_fht(n, M, L, ordering=None)
    def Ab(b):
        return Ax(b).reshape(-1, 1) / np.sqrt(n)
    def Az(z):
        return Ay(z).reshape(-1, 1) / np.sqrt(n)
    return Ab, Az, ordering
# return the ordering
# have second function which does the same but uses the previous ordering. And on fewer sections. 
# Check it on some betas. Use same beta and check it's the same result. Set last sections of beta to zero 
# For the second function just use the first L unprotected rows of the ordering. 


#SPARC dictionary
#Returns 2 functions Ab and Az which compute Aβ and A^Tz respectively.
#Takes in the value of ordering so a smaller design matrix can be generated to run on 
#just the unprotected sections for the second round of amp decoding. 
def sparc_transforms_shorter(L, M, n, ordering):
    # just use the first L unprotected rows of ordering
    Ax, Ay, _ = block_sub_fht(n, M, L, ordering=ordering[:L,:])
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
        #if τ == last_τ:
        # can change params for isclose if need be. If stopping too early.
        if np.isclose(τ, last_τ):
            return β
        last_τ = τ
        
        # added astype to increase precision to avoid divide by zero in LLR
        s = (β + Az(z)).astype(np.longdouble)
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
	def __init__(self, L, M, sigma, p, r, t):
		self.L = L
		self.M = M
		self.sigma = sigma
		self.p = p
		self.r = r
		#self.r_pa = r_pa
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
    #note L may not equal the total number of sections for the sparc code. 
    #normally just interested in the number of sections covered by the LDPC
    
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


# function to convert bitwise probabilities to sectionwise probabilities
# v is the bitwise probabilities
# L is the number of sections in v
# M is the number of bits in v
def bp2sp(v, L, M: "must be power of 2"):
    logm = int(np.log2(M))
    # sectionwise probabilities
    sp = np.zeros(L*M)

    for l in range(L):
        # bitwise probabilities for section l
        bp = v[l*logm:(l+1)*logm]

        for m in range(M):
            # The bit values that correspond to the mth entry of this section being non-zero
            bits = bitarray(bin(m)[2:].zfill(logm))
            # bp and (1-bp) to the power of bits ensures that only the probability of each bit 
            # being the value we would need it to be for this column to be non-zero is included
            # in the probability. 
            # So we are multiplying the probability of each bit being how we need it to make the
            # current column non-zero.
            print("bp", bp)
            print("bits", bits)
            a = bp**bits
            bits.invert()
            b = (1-bp)**(bits)
            sp[l*M+m] = np.prod(a * b)

    return sp

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

# note removed r_pa from function as never use it
def amp_ldpc_sim(sparcparams: SPARCParams, ldpcparams: LDPCParams = None):
    #Get the SPARC parameters from the struct
    L = sparcparams.L
    M = sparcparams.M
    P = sparcparams.p
    sigma = sparcparams.sigma
    r_sparc = sparcparams.r
    #r_pa_sparc = sparcparams.r_pa
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
    unprotected_bits = np.random.randint(0, 2, int(total_bits-nl)).tolist()

    # concatenate the ldpc and unprotected bits      
    sparc_bits = unprotected_bits+ldpc_bits

    assert len(sparc_bits)==total_bits

    # convert the bits to indices for encoding
    sparc_indices = bits2indices(sparc_bits, M)

    assert len(sparc_indices)==L
       
    # Generate the SPARC transform functions A.beta and A'.z
    Ab, Az, ordering = sparc_transforms(L, M, n)
    
    # Generate our transmitted signal X
    β_0 = np.zeros((L*M, 1))
    for l in range(L):
        β_0[l*M + sparc_indices[l]] = np.sqrt(n * Pl[l])
    x = Ab(β_0)
    
    # check that the power has been allocated uniformly. This should be approx equal to one when divided by the snr. 
    #print("Average power/snr is ", np.mean(x**2)/snr)

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
    #correct_amp = np.sum(np.array(rx_message) == np.array(sparc_indices)) / L
    #print("Fraction of sections decoded correctly with amp: ", correct_amp)
    
    # Compute BER for just sparc
    ber_amp = sum(bin(a^b).count('1') for (a, b) in zip(sparc_indices, rx_message))/total_bits


    if ldpcparams==None:
        ber_ldpc = None
        ber_ldpc_amp = None
    else:
        # Get the sectionwise posterior probabilities by dividing β by the power in each section. 
        # This converts each section in β to a valid probability distribution.
        sectionwise_posterior = β/np.sqrt(n*P/L)

        ldpc_sections = int(nl/logm)
        # convert sectionwise to bitwise posterior probabilities
        # only need the sections corresponding to ldpc code 
        bitwise_posterior = sp2bp(sectionwise_posterior[(L-ldpc_sections)*M:], ldpc_sections, M) 

        np.clip(bitwise_posterior, 0.001, 1-0.001, out=bitwise_posterior)
        # computer the log likelihood ratio for decoding
        LLR = np.log(1-bitwise_posterior)- np.log(bitwise_posterior)

        (app, it) = ldpc_code.decode(LLR)

        test_output = (bitwise_posterior>0.5)
        #print("Incorrect bits before ldpc decoding: ", np.sum(test_output!=ldpc_bits))

        # check that the inequality is the right way round in below
        v_output = (app<0.0)
        #print("Incorrect bits in ldpc decoding: ", np.sum(v_output!=ldpc_bits))

        # convert the bits to indices
        beta_output_ldpc = bits2indices(v_output, M)

        #### Compute performance with just amp and then ldpc
        beta_output = rx_message
        beta_output[L-ldpc_sections:] = beta_output_ldpc
        # compute fraction of sections decoded correctly with AMP followed by SPARC decoding
        #correct_ldpc = np.sum(np.array(beta_output) == np.array(sparc_indices))/L
        #print("Fractions of sections decoded correctly with amp and ldpc: ", correct_ldpc)

        # compute BER for amp and ldpc
        ber_ldpc = sum(bin(a^b).count('1') for (a, b) in zip(sparc_indices, beta_output))/total_bits

        # only perform final round of AMP decoding if not all sections of beta are covered with LDPC
        if (L-ldpc_sections==0):
            ber_ldpc_amp = None
        else:    
            # generate beta_ldpc which has the first L unprotected sections set to zero
            # and the L ldpc sections will have exactly one non-zero entry per section. 
            beta_ldpc = np.zeros((L*M, 1))
            # using variable i because beta_output_L only has ldpc_sections number of entries. 
            i=0
            for l in range(L-ldpc_sections,L):
                beta_ldpc[(l)*M + beta_output_ldpc[i]] = np.sqrt(n * Pl[l])
                i+=1
            x_ldpc = Ab(beta_ldpc)

            # calculate a new value of the channel output which doesn't contain the contribution from 
            # the ldpc sections.
            y_new = y - x_ldpc

            # Rerun amp decoding on the new input y_new over the unprotected sections.
            L_unprotected = L-ldpc_sections
            Ab_new, Az_new = sparc_transforms_shorter(L_unprotected, M, n, ordering)

            beta_new = amp(y_new, sigma, Pl[:L_unprotected], L_unprotected, M, T, Ab_new, Az_new).reshape(-1)

            # Convert decoded beta_new back to a message
            # beta_new gives you the first L-ldpc_sections sections of the final received message. 
            rx_message_final = []
            for l in range(L_unprotected):
                idx = np.argmax(beta_new[l*M:(l+1)*M])
                rx_message_final.append(idx)


            #beta_output = np.zeros(L)
            beta_output[:L_unprotected] = rx_message_final
            # don't need to do line below as they're already equal
            #beta_output[L-ldpc_sections:] = beta_output_ldpc
            # compute fraction of sections decoded correctly with AMP followed by SPARC decoding followed by amp.
            #correct_ldpc = np.sum(np.array(beta_output) == np.array(sparc_indices))/L
            #print("Fractions of sections decoded correctly with amp and ldpc and amp: ", correct_ldpc)

                
            # compute BER for amp and ldpc and amp again.
            ber_ldpc_amp = sum(bin(a^b).count('1') for (a, b) in zip(sparc_indices, beta_output))/total_bits
        
            

    # overall rate of the code
    R = (L*logm - (nl-kl))/n

    #Compute Eb/N0
    #EbN0 = 1/(2*R) * (P/sigma**2)
   
    return ber_amp, ber_ldpc, ber_ldpc_amp, R




if __name__ == "__main__":
    # get the time so you can calculate the wall clock time of the process
    t0 = time.time()
    '''
    ##########################################
    # testing bp2sp function
    L=2
    M=4
    # sectionwise probabilities
    beta = [0.8, 0.1, 0.05, 0.05, 0.0, 0.05, 0.1, 0.85]
    print("original sectionwise: ", beta)
    # convert this to bitwise probs
    v = sp2bp(beta, L, M)
    print("bitwise: ", v)
    sp = bp2sp(v, L, M)
    print("sectionwise: ", sp)
    '''
    '''

    ########################################################
    # Keep rate fixed and vary the EbN0
    # plot of performance of sparc with outer code and amp only, after LDPC, after final AMP
    # SPARC without outer code. 
    # sparc with outer code and amp only, will have a higher sparc rate than the overall and won't be getting the benefits from ldpc so will therefore perform worse
    # sparc without outer code with have the same overall rate as the sparc with the LDPC applied so is a better comparison 
    
    
    #Compute Eb/N0
    #EbN0 = 1/(2*R) * (P/sigma**2)
    # Sparc parameters
    L = 768
    M = 512
    logm = np.log2(M)
    p=1.5
    r_sparc = 1
    T = 64

    logm = np.log2(M)

    # LDPC parameters
    standard = '802.16'
    r_ldpc = '5/6'
    # number of ldpc sections, must be divisible by 8 to ensure nl is divisible by 24.
    # covering roughly 73% of sections like in Adams paper
    sec = 768
    # number of ldpc bits, must be divisible by 24
    nl = logm * sec
    z = int(nl/24)
    ldpcparams = LDPCParams(standard, r_ldpc, z)    

    n = L*logm/r_sparc
    R = (L*logm-nl*(1-5/6))/n
    print('Overall rate is: ', R)
    # need to change the 5/6 in here if I change r_ldpc.            
    # the Eb/N0 at capacity - I don't think this worked.
    #snrc = (2**(2*R) - 1)
    #EbN0c = 1/(2*R) * snrc
    #EbN0c_dB = 20*log10(EbN0c)
    # capacity of sigma



    repeats = 1
    datapoints = 8
    SIGMA = linspace(1, 0.3, datapoints)
    BER_amp = np.zeros(datapoints)
    BER_ldpc = np.zeros(datapoints)
    BER_ldpc_amp = np.zeros(datapoints)
    BER_sparc = np.zeros(datapoints)

    i=0
    for sigma in SIGMA:
        sparcparams_outer = SPARCParams(L, M, sigma, p, r_sparc, T)
        BER_ldpc_rep = np.zeros(repeats)

        for j in range(repeats):
            (_, BER_ldpc_rep[j], _, _) = amp_ldpc_sim(sparcparams_outer, ldpcparams)

        BER_ldpc[i] = np.sum(BER_ldpc_rep)/repeats
        i+=1
    #snrdB = 20*np.log10(P/sigma**2)
    EbN0 = 1/(2*R) * (p/SIGMA**2)
    print(EbN0)
    EbN0_dB = 20*log10(EbN0)
    # open file you want to write CSV output to. 'a' means its in append mode. Switching this to 'w' will make it overwrite the file.
    myFile = open('EbN0_dBVsBER_amp_ldpc_1.csv', 'a')
    with myFile:
        myFields = ['EbN0_dB', 'BER_ldpc']
        writer = csv.DictWriter(myFile, fieldnames=myFields)
        writer.writeheader()
        for k in range(datapoints):
            writer.writerow({'EbN0_dB' : EbN0_dB[k], 'BER_ldpc' : BER_ldpc[k]})
    

    fig, ax = plt.subplots()
    ax.set_yscale('log', basey=10)
    ax.plot(EbN0_dB, BER_ldpc, 'k:', label="SPARC with outer code: after LDPC")
    #plt.axvline(x=EbN0c_dB, color='r', linestyle='-', label='Shannon limit')
    plt.xlabel('EbN0 in dB') # at some point need to work out how to write this so it outputs properly
    plt.ylabel('BER')
    #plt.title("All sections of code covered. R_sparc=1, R_ldpc=5/6")
    plt.legend()
    print("Wall clock time elapsed: ", time.time()-t0)
    plt.show()
    plt.savefig('EbN0VsBER_amp_ldpc_1.png')

    #########################################################
    '''
    ''' 
    ########################################
    Delete this code soon. Just saving incase I want to recheck it when I'm less sleepy. 
    # testing to see if new, smaller design matrix and Az work
    L = 1024
    L_ldpc = 256
    M = 512
    logm = int(np.log2(M))
    r_sparc = 1
    n = int(L*np.log2(M) / r_sparc)
    L_unprotected = L-L_ldpc
    P=1

    # Generate the remaining bits required
    bits = np.random.randint(0, 2, L*logm).tolist()

    # convert the bits to indices for encoding
    beta = bits2indices(bits, M)

    Ab, Az, ordering = sparc_transforms(L, M, n)

    # want to keep the value of n the same as this determines the number of rows and we want to keep the number of rows constant
    # It is just giving us an effectively lower rate
    Ab_u, Az_u = sparc_transforms_shorter(L_unprotected, M, n, ordering)
    Pl = P/L * np.ones(L)
     # Generate our transmitted signal X
    β_0 = np.zeros((L*M, 1))
    for l in range(L):
        β_0[l*M + beta[l]] = np.sqrt(n * Pl[l])


    β_0[(L-L_ldpc)*M:] = np.zeros(L_ldpc*M).reshape(-1,1)
    beta_2 = β_0[:(L-L_ldpc)*M]
    

    assert(sum(β_0)==sum(beta_2))  
    x = Ab(β_0)
    x2 = Ab_u(beta_2)

    size = x2.size
    print(x[:size]-x2)
    assert((x[:size]==x2).all())
    assert(sum(x[size:])==0)

    a = Az(x)
    b = Az_u(x2)
    size2 = b.size
    print(sum(a[size2:]))
    print(a[:size2]-b)
    assert((a[:size2]==b).all())
    # this final assertion doesn't matter. We don't care what Az() does to the final values.
    # as we want to essentially force these values to zero. As long as the first lot of values are the same
    assert(sum(a[size2:])==0)
   


    # get the time so you can calculate the wall clock time of the process
    t0 = time.time()
    '''
    '''c = np.array([0.2, 0.3, 0.1, 0.4, 0.6, 0.4, 0, 0])
    print(sp2bp(c, 2, 4))

    ##############################################################
    # Number of sections covered by ldpc code vs. BER.    
    
    '''
    '''M = 512
    L = 512
    sigma = 1
    P = 3
    r_sparc = 0.75
    # removed r_pa from function as never us it
    #r_pa_sparc = 1
    T = 64
    sparcparams = SPARCParams(L, M, sigma, P, r_sparc, T)
    
    # The ldpc sections must be less than L and must be divisible by 8 to give an integer value for z
    ldpc_sections = 400
    # nl is the number of bits that will be LDPC bits. 
    # This is given by the no. bits per sections * the no. sections
    # Note the no. sections must be at least 8 to give something divisible by 24
    nl = np.log2(512) * ldpc_sections
    # Can get from nl to z. z determines the size of the protograph matrix
    z = int(nl/24)
    standard = '802.16'
    rate = '5/6'
    ldpcparams = LDPCParams(standard, rate, z)

    (ber_amp, ber_ldpc, R) = amp_ldpc_sim(sparcparams)#, ldpcparams)
    print("ber_amp: ", ber_amp)
    #print("ber_ldpc: ", ber_ldpc)
    print("R: ", R)

    datapoints = 9
    # This should give increments of 64 each time
    # Ensures number of sections is divisible by 8. 
    ldpc_sections = linspace(0, 512, datapoints)

    ber_ldpc = np.zeros(datapoints)
    R_overall = np.zeros(datapoints)
    i = 0
    repeats = 10
    M = 512
    L = 1024
    logm = int(np.log2(M))
    # Overall rate R
    R = 0.75
    for sec in ldpc_sections:
        #number of LDPC bits. Must be divisible by 8.
        nl = logm * sec
        print(nl)
        # Can get from nl to z. z determines the size of the protograph matrix
        z = int(nl/24)
        standard = '802.16'
        # if change this, change formula for kl
        rate = '2/3'
        
        ldpcparams = LDPCParams(standard, rate, z)
        # note 2/3 is the rate of the ldpc. Hard coded this to avoid rounding errors in storing this
        kl = nl * 2/3
        n = (L*logm - (nl-kl))/R
        # overall rate of the ldpc/ sparc code. 
        r_sparc = (L*logm)/n


        sigma = 1
        P = 2.3
        T = 64
        # sparc params for the ldpc/ sparc code
        sparcparams_ldpc = SPARCParams(L, M, sigma, P, r_sparc, T)

        ber_ldpc_rep = np.zeros(repeats)
        R_overall_rep = np.zeros(repeats)
        for j in range(repeats):
            if z>0:
                (_, ber_ldpc_rep[j], R_overall_rep[j]) = amp_ldpc_sim(sparcparams_ldpc, ldpcparams)
            else:
                # if z=0 then just doing sparc with no ldpc
                (ber_ldpc_rep[j], _, R_overall_rep[j]) = amp_ldpc_sim(sparcparams_ldpc)
        ber_ldpc[i] = np.sum(ber_ldpc_rep)/repeats
        R_overall[i] = np.sum(R_overall_rep)/repeats
        i += 1
    print(R_overall)
    plot(ldpc_sections, ber_ldpc)
    xlabel('Sections covered by ldpc')
    ylabel('BER')
    show()'''
    '''
    ################################################################
    # Keep the sparc rate constant and then the snr we choose will give the same
    # BER each time. But vary the number of sections covered. This will vary the overall rate
    # Compare to a sparc code with the overall rate and see which performs better. 
    datapoints = 17
    # This should give increments of 64 each time
    # Ensures number of sections is divisible by 8. 
    ldpc_sections = linspace(0, 1024, datapoints)
    print(ldpc_sections)

    ber_ldpc = np.zeros(datapoints)
    ber_noldpc = np.zeros(datapoints)
    R_overall = np.zeros(datapoints)
    i = 0
    repeats = 1
    M = 512
    L = 1024
    logm = int(np.log2(M))
    # sparc rate
    r_sparc = 0.75
    for sec in ldpc_sections:
        #number of LDPC bits. Must be divisible by 8.
        nl = logm * sec
        print(nl)
        # Can get from nl to z. z determines the size of the protograph matrix
        z = int(nl/24)
        standard = '802.16'
        # if change this, change formula for kl
        rate = '5/6'
        
        ldpcparams = LDPCParams(standard, rate, z)

        sigma = 1
        P = 2.4
        T = 64
        # sparc params for the ldpc/ sparc code
        sparcparams_ldpc = SPARCParams(L, M, sigma, P, r_sparc, T)
        ber_ldpc_rep = np.zeros(repeats)
        ber_noldpc_rep = np.zeros(repeats)
        R_overall_rep = np.zeros(repeats)
        for j in range(repeats):
            if z>0:
                (_, ber_ldpc_rep[j], R_overall_rep[j]) = amp_ldpc_sim(sparcparams_ldpc, ldpcparams)
                # calculate the BER for a sparc with the same overall rate but no ldpc code.
                sparcparams_noldpc = SPARCParams(L, M, sigma, P, R_overall_rep[j], T)
                (ber_noldpc_rep[j],_,_) = amp_ldpc_sim(sparcparams_noldpc)
            else:
                # if z=0 then just doing sparc with no ldpc
                (ber_ldpc_rep[j], _, R_overall_rep[j]) = amp_ldpc_sim(sparcparams_ldpc)
                ber_noldpc_rep[j] = ber_ldpc_rep[j]
        ber_ldpc[i] = np.sum(ber_ldpc_rep)/repeats
        R_overall[i] = np.sum(R_overall_rep)/repeats
        ber_noldpc[i] = np.sum(ber_noldpc_rep)/repeats
        i += 1

    fig, ax = plt.subplots()
    ax.set_yscale('log', basey=10)
    ax.plot(R_overall, ber_ldpc, 'k--', label = 'With ldpc')
    ax.plot(R_overall, ber_noldpc, 'k:', label = 'Without ldpc')
    plt.xlabel('Overall Rate')
    plt.ylabel('BER')
    plt.legend()
    plt.title('Trade off with using lower overall rate to cover more sections with ldpc or just having lower sparc rate.')
    plt.show()
    print("Wall clock time elapsed: ", time.time()-t0)
    '''
    
    '''
    ########################################################
    # keep snr fixed and vary the rate
    # plot of performance of sparc with outer code and amp only, after LDPC, after final AMP
    # SPARC without outer code. 
    # sparc with outer code and amp only, will have a higher sparc rate than the overall and won't be getting the benefits from ldpc so will therefore perform worse
    # sparc without outer code with have the same overall rate as the sparc with the LDPC applied so is a better comparison 

    # Sparc parameters
    L = 1024
    M = 512
    sigma = 1
    # pick one of the 3 powers below 
    P = 2.4
    T = 64
    # the rate at capacity for P and sigma
    capacity = 1/2*np.log2(1 + P/sigma**2)
    logm = np.log2(M)

    # LDPC parameters
    standard = '802.16'
    r_ldpc = '5/6'
    # number of ldpc sections, must be divisible by 8 to ensure nl is divisible by 24.
    # covering roughly 73% of sections like in Adams paper
    sec = 768
    # number of ldpc bits, must be divisible by 24
    nl = logm * sec
    z = int(nl/24)
    ldpcparams = LDPCParams(standard, r_ldpc, z)
        

    repeats = 100
    datapoints =7
    R = linspace(0.3, 0.6, datapoints)
    BER_amp = np.zeros(datapoints)
    BER_ldpc = np.zeros(datapoints)
    BER_ldpc_amp = np.zeros(datapoints)
    BER_sparc = np.zeros(datapoints)

    i=0
    for r in R:
        # need to change the 5/6 in here if I change r_ldpc.
        n = (L*logm-nl*(1-5/6))/r
        r_sparc = L*logm/n 
        sparcparams_outer = SPARCParams(L, M, sigma, P, r_sparc, T)
        sparcparams = SPARCParams(L, M, sigma, P, r, T)
        BER_amp_rep = np.zeros(repeats)
        BER_ldpc_rep = np.zeros(repeats)
        BER_ldpc_amp_rep = np.zeros(repeats)
        BER_sparc_rep = np.zeros(repeats)

        for j in range(repeats):
            (BER_amp_rep[j], BER_ldpc_rep[j], BER_ldpc_amp_rep[j], _) = amp_ldpc_sim(sparcparams_outer, ldpcparams)
            (BER_sparc_rep[j],_,_,_) = amp_ldpc_sim(sparcparams, ldpcparams)

        BER_amp[i] = np.sum(BER_amp_rep)/repeats
        BER_ldpc[i] = np.sum(BER_ldpc_rep)/repeats
        BER_ldpc_amp[i] = np.sum(BER_ldpc_amp_rep)/repeats
        BER_sparc[i] = np.sum(BER_sparc_rep)/repeats
        i+=1

    # open file you want to write CSV output to. 'a' means its in append mode. Switching this to 'w' will make it overwrite the file.
    myFile = open('RVsBER_amp_ldpc_3.csv', 'a')
    with myFile:
        myFields = ['R', 'BER_amp', 'BER_ldpc', 'BER_ldpc_amp', 'BER_sparc']
        writer = csv.DictWriter(myFile, fieldnames=myFields)
        writer.writeheader()
        for k in range(len(R)):
            writer.writerow({'R' : R[k], 'BER_amp' : BER_amp[k], 'BER_ldpc' : BER_ldpc[k], 'BER_ldpc_amp' : BER_ldpc_amp[k], 'BER_sparc' : BER_sparc[k]})
    

    fig, ax = plt.subplots()
    ax.set_yscale('log', basey=10)
    ax.plot(R, BER_amp, 'k--', label = 'SPARC with outer code: AMP only')
    ax.plot(R, BER_ldpc, 'k:', label = 'SPARC with outer code: after LDPC')
    ax.plot(R, BER_ldpc_amp, 'k-', label = 'SPARC with outer code: after final AMP')
    ax.plot(R, BER_sparc, 'k-.', label = 'SPARC with no outer code')
    plt.axvline(x=capacity, color='r', linestyle='-', label='Shannon limit')
    plt.xlabel('Overall Rate')
    plt.ylabel('BER')
    plt.legend()
    print("Wall clock time elapsed: ", time.time()-t0)
    plt.savefig('RVsBER_amp_ldpc_3.png')
    '''

    
    ########################################################
    # Keep rate fixed and vary the snr
    # plot of performance of sparc with outer code and amp only, after LDPC, after final AMP
    # SPARC without outer code. 
    # sparc with outer code and amp only, will have a higher sparc rate than the overall and won't be getting the benefits from ldpc so will therefore perform worse
    # sparc without outer code with have the same overall rate as the sparc with the LDPC applied so is a better comparison 
    
    
    #Compute Eb/N0
    #EbN0 = 1/(2*R) * (P/sigma**2)
    # Sparc parameters
    L = 768
    M = 512
    logm = np.log2(M)
    sigma = 1
    r_sparc = 0.94
    T = 64

    logm = np.log2(M)

    # LDPC parameters
    standard = '802.16'
    r_ldpc = '5/6'
    # number of ldpc sections, must be divisible by 8 to ensure nl is divisible by 24.
    # covering roughly 73% of sections like in Adams paper
    sec = 568
    # number of ldpc bits, must be divisible by 24
    nl = logm * sec
    z = int(nl/24)
    ldpcparams = LDPCParams(standard, r_ldpc, z)    

    n = L*logm/r_sparc
    R = (L*logm-nl*(1-5/6))/n
    print('Overall rate is: ', R)
    # need to change the 5/6 in here if I change r_ldpc.            
    # the Eb/N0 at capacity 
    snrcdB = 20*np.log10(2**(2*R) - 1)
    #EbN0c = 1/(2*R) * snrc

    repeats = 1
    datapoints = 8
    P = linspace(1.5, 5.5, datapoints)
    BER_amp = np.zeros(datapoints)
    BER_ldpc = np.zeros(datapoints)
    BER_ldpc_amp = np.zeros(datapoints)
    BER_sparc = np.zeros(datapoints)

    i=0
    for p in P:
        sparcparams_outer = SPARCParams(L, M, sigma, p, r_sparc, T)
        sparcparams = SPARCParams(L, M, sigma, p, R, T)
        BER_amp_rep = np.zeros(repeats)
        BER_ldpc_rep = np.zeros(repeats)
        BER_ldpc_amp_rep = np.zeros(repeats)
        BER_sparc_rep = np.zeros(repeats)

        for j in range(repeats):
            (BER_amp_rep[j], BER_ldpc_rep[j], BER_ldpc_amp_rep[j], _) = amp_ldpc_sim(sparcparams_outer, ldpcparams)
            (BER_sparc_rep[j],_,_,_) = amp_ldpc_sim(sparcparams, ldpcparams)

        BER_amp[i] = np.sum(BER_amp_rep)/repeats
        BER_ldpc[i] = np.sum(BER_ldpc_rep)/repeats
        BER_ldpc_amp[i] = np.sum(BER_ldpc_amp_rep)/repeats
        BER_sparc[i] = np.sum(BER_sparc_rep)/repeats
        i+=1
    snrdB = 20*np.log10(P/sigma**2)
    #EbN0 = 1/(2*R) * (P/sigma**2)
    # open file you want to write CSV output to. 'a' means its in append mode. Switching this to 'w' will make it overwrite the file.
    myFile = open('SNR_dBVsBER_amp_ldpc_2.csv', 'a')
    with myFile:
        myFields = ['SNR_dB', 'BER_amp', 'BER_ldpc', 'BER_ldpc_amp', 'BER_sparc']
        writer = csv.DictWriter(myFile, fieldnames=myFields)
        writer.writeheader()
        for k in range(datapoints):
            writer.writerow({'SNR_dB' : snrdB[k], 'BER_amp' : BER_amp[k], 'BER_ldpc' : BER_ldpc[k], 'BER_ldpc_amp' : BER_ldpc_amp[k], 'BER_sparc' : BER_sparc[k]})
    

    fig, ax = plt.subplots()
    ax.set_yscale('log', basey=10)
    ax.plot(snrdB, BER_amp, 'k--', label = 'SPARC with outer code: AMP only')
    ax.plot(snrdB, BER_ldpc, 'k:', label = 'SPARC with outer code: after LDPC')
    ax.plot(snrdB, BER_ldpc_amp, 'k-', label = 'SPARC with outer code: after final AMP')
    ax.plot(snrdB, BER_sparc, 'k-.', color='g', label = 'SPARC with no outer code')
    plt.axvline(x=snrcdB, color='r', linestyle='-', label='Shannon limit')
    plt.xlabel('SNR') # at some point need to work out how to write this so it outputs properly
    plt.ylabel('BER')
    plt.legend()
    print("Wall clock time elapsed: ", time.time()-t0)
    plt.savefig('SNR_dBVsBER_amp_ldpc_2.png')
    
    '''
    #########################################################
    # plot of performance of sparc code on it's own with a decreasing rate. 
    # Sparc parameters
    L = 1024
    M = 512
    sigma = 1
    # 3 different values of power to test. 
    P1 = 2.0
    P2 = 2.2 
    P3 = 2.4
    T = 64

    repeats = 1
    datapoints = 6
    R = linspace(0.5, 0.75, datapoints)
    BER1 = np.zeros(datapoints)
    BER2 = np.zeros(datapoints)
    BER3 = np.zeros(datapoints)
    i=0
    for r_sparc in R:
        sparc_params1 = SPARCParams(L, M, sigma, P1, r_sparc, T)
        sparc_params2 = SPARCParams(L, M, sigma, P2, r_sparc, T)
        sparc_params3 = SPARCParams(L, M, sigma, P3, r_sparc, T)
        BER1_rep = np.zeros(repeats)
        BER2_rep = np.zeros(repeats)
        BER3_rep = np.zeros(repeats)
        for j in range(repeats):
            (BER1_rep[j],_,_,_) = amp_ldpc_sim(sparc_params1)
            (BER2_rep[j],_,_,_) = amp_ldpc_sim(sparc_params2)
            (BER3_rep[j],_,_,_) = amp_ldpc_sim(sparc_params3)
        BER1[i] = np.sum(BER1_rep)/repeats
        BER2[i] = np.sum(BER2_rep)/repeats
        BER3[i] = np.sum(BER3_rep)/repeats
        i+=1
    fig, ax = plt.subplots()
    ax.set_yscale('log', basey=10)
    ax.plot(R, BER1, 'k--', label = 'Power=2.0')
    ax.plot(R, BER2, 'k:', label = 'Power=2.2')
    ax.plot(R, BER3, 'k-.', label = 'Power=2.4')
    plt.xlabel('Overall Rate')
    plt.ylabel('BER')
    plt.legend()
    plt.title("BER of the SPARC code as rate varies")
    print("Wall clock time elapsed: ", time.time()-t0)
    plt.savefig('sparc_RVsBER_1.png')
    '''
    '''
    #################################################################
    #Code to plot the performance of the sparc code on it's own for 
    # a fixed rate and varying snr. 
    L = 1024
    M = 512
    # fix rate
    R = 0.75

    σ_n = 1.0
    # shannon limit on snr. i.e. snr at capacity, where R is capacity.
    snrc = 2**(2*R)-1
    
    # Unsure what values to choose for T and R_PA. Just copied previous example.
    T = 64
    # store the error rate for each snr
    error_rate = []
    repeats = 1
    snr_array = np.linspace(snrc, snrc+1, 11)
    print(snr_array)
    for snr in snr_array:
        # calculate the new value of P to adjust the snr
        P = snr * σ_n**2
        error = []
        sparcparams = SPARCParams(L, M, σ_n, P, R, T)
        for i in range(repeats):
            ber = amp_ldpc_sim(sparcparams)[0]
            error.append(ber)
        aver_error = np.sum(error)
        error_rate.append(aver_error/repeats)
        
    plot(snr_array, error_rate)
    xlabel('snr')
    title("Parameters: L=1024, M=512, σ_n=1, R=0.75, T=64")
    ylabel('BER')
    show()'''


