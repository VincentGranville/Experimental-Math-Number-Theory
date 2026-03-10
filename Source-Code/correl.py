# Compute binary digits of X, p*X, q*X backwards (assuming X is random)
# Only digits after the decimal point (on the right) are computed
# Compute correlations between digits of p*X and q*X
# Include carry-over when performing grammar school multiplication

import numpy as np
import gmpy2

kmax = 10000000
ctx = gmpy2.get_context()  
ctx.precision = kmax 
ndigits = ctx.precision  
z = gmpy2.sqrt(2)

# main parameters
seed = 195
np.random.seed(seed)
# p, q odd integers, coprime
p  = 3   
q  = 5  

def gmpy2_correl(z, p, q):

    # correl b/w binary digits of z and pz/q (needs p < q)
    zstri = gmpy2.digits(z, 2)[0]   # get binary digits of z as a string
    zoff = gmpy2.digits(z, 2)[1]

    w = gmpy2.mpfr(z*p)/gmpy2.mpz(q) 
    woff = gmpy2.digits(w, 2)[1]
    w_offset = '0' * (zoff - woff)  # works only if p < q
    wstri = w_offset + gmpy2.digits(w, 2)[0]

    prod = 0 
    for k in range(kmax):
        d1 = int(zstri[k])
        d2 = int(wstri[k])
        prod +=  d1*d2
        correl = 4*prod/(k+1) - 1
        if k % 100000 == 0 and k > 100:
            checksum = correl * p * q  # should be close to 1
            print("gmpy2> k: %7d correl: %9.7f check: %9.7f" %(k, correl, checksum))
    return(correl)


def vg_correl(z, p, q):

    # correl b/w binary digits of pz and qz
    mode = 'constant'  # options: 'random', 'constant'

    # local variables
    zstri = gmpy2.digits(z, 2)[0]   # get binary digits of z as a string
    X, pX, qX = 0, 0, 0
    d1, d2, e1, e2 = 0, 0, 0, 0
    prod, count = 0, 0 
    sum1 = 0
    sum2 = 0

    # loop over digits in reverse order
    for k in range(kmax):

        # b is a digit of X
        if mode == 'random':
            b = np.random.randint(0, 2)  
        else:
            b = int(zstri[kmax-k-1])
        X = b + X/2  

        c1 = p*b
        old_d1 = d1
        old_e1 = e1 
        d1 = (c1 + old_e1//2) %2  # digit of pX
        e1 = (old_e1//2) + c1 - d1
        pX = d1 + pX/2

        c2 = q*b
        old_d2 = d2
        old_e2 = e2 
        d2 = (c2 + old_e2//2) %2  #digit of qX
        e2 = (old_e2//2) + c2 - d2
        qX = d2 + qX/2

        prod  += d1*d2
        count += 1 
        sum1 += d1
        sum2 += d2
        mean1 = sum1/count
        mean2 = sum2/count
        std1 = (mean1 * (1 - mean1))**0.5 
        std2 = (mean2 * (1 - mean2))**0.5 
        covar = prod/count - mean1*mean2
        if count > 100:
            correl = covar/(std1*std2)
        else:
            correl = 0 
        #correl = 4*prod/count - 1 

        if k% 100000 == 0:  
            checksum = p*q*correl # should be close to 1
            print("vg>k = %7d, correl = %9.6f checksum = %9.6f" % (k, correl, checksum))  

    print("\np = %3d, q = %3d" %(p, q))
    print("X = %12.9f, pX  = %12.9f, qX  = %12.9f" % (X, pX, qX))
    print("X = %12.9f, p*X = %12.9f, q*X = %12.9f" % (X, p*X, q*X))    
    print("Correl = %7.4f, 1/(p*q) = %7.4f" % (correl, 1/(p*q))) 
    return(correl)

#--- Main

correl1 = gmpy2_correl(z, p, q)
correl2 = vg_correl(z, p, q)

