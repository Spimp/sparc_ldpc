# Code removed from the main file in sparc_ldpc.py. This contains a lot of the code
# for plotting graphs. Trying to transfer this over to functions so the main function is kept
# clearer. 

    ######################################
    # Plot the parameterised power allocation
    datapoints = 5
    P=2.5
    L=512
    R=1
    SNR_dB = np.linspace(10, 14, datapoints)
    SNR = 10**(SNR_dB/20)
    pa = np.zeros((datapoints,L))
    i=0
    for snr in SNR:
        # The channel capacity of the AWGN channel 
        #snr = P / sigma**2
        C = 0.5 * np.log2(1 + snr)
        # A good starting point for the parameters is:
        a=R/C
        f=R/C
        pa[i,:] = pa_parameterised(L, C, P, a, f)
        i=i+1
    fig, ax = plt.subplots()
    ax.plot(pa[0,:], 'b--', label='$SNR_{dB}$='+str(SNR_dB[0])) 
    ax.plot(pa[1,:], 'k--', label='$SNR_{dB}$='+str(SNR_dB[1]))
    ax.plot(pa[2,:], 'm--', label='$SNR_{dB}$='+str(SNR_dB[2]))
    ax.plot(pa[3,:], 'c--', label='$SNR_{dB}$='+str(SNR_dB[3]))
    ax.plot(pa[4,:], 'r--', label='$SNR_{dB}$='+str(SNR_dB[4]))

    plt.xlabel('L')
    plt.ylabel('Power Allocation')
    plt.legend(loc=6, prop={'size': 7})
    plt.title("The power allocation for a rate 1 SPARC with different SNRs")
    plt.show()

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
    # Keep rate fixed and vary the EbN0 by varying sigma to give a waterfall 
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
    p=1.8
    r_sparc = 1
    T = 64

    logm = np.log2(M)

    # LDPC parameters
    standard = '802.16'
    r_ldpc = '5/6'
    # number of ldpc sections, must be divisible by 8 to ensure nl is divisible by 24.
    # covering all sections to give overall rate of 5/6
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
    # capacity of sigma. Calculated using J(2/sigma^2)=5/6
    # The approximations in the appendix for paper "Design of LDPC codes..." was used
    sigma2c = 2/(-0.706692*log(0.386013*(1-5/6))+1.75017*(5/6))
    EbN0c = 1/(2*R) * (p/sigma2c)
    EbN0c_dB = 10*log10(EbN0c)


    repeats = 100
    datapoints = 7
    SIGMA = linspace(1, 0.4, datapoints)
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
    ax.plot(EbN0_dB, BER_ldpc, 'k:', label = 'SPARC with outer code: after LDPC')
    plt.axvline(x=EbN0c_dB, color='r', linestyle='-', label='Shannon limit')
    plt.xlabel('EbN0_dB') # at some point need to work out how to write this so it outputs properly
    plt.ylabel('BER')
    plt.tight_layout()
    plt.legend()
    print("Wall clock time elapsed: ", time.time()-t0)
    plt.savefig('EbN0_dBVsBER_amp_ldpc_1.png')

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
    myFile = open('SNR_dBVsBER_amp_ldpc_test.csv', 'a')
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
    plt.savefig('SNR_dBVsBER_amp_ldpc_test.png')
    '''
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


