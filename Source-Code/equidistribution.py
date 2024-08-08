import numpy as np
import random
seed = 1
random.seed(seed) # try 0, 1, 2

def int_to_binstring(x, n):
    # convert integer x into n-bit string
    str = bin(x).replace("0b",'')
    while len(str) < n:
        str = '0' + str
    return(str)

def hash_update(key, hash):
    if key in hash:
        hash[key] += 1
    else:
        hash[key] = 1
    return(hash)

def summary(hash, m, N):
    # summary stats for congruences modulo m
    mean = 0
    std = 0
    for j in range(m):
        # loop on residue classes modulo m
        if (m, j) in hash:
            count = hash[(m, j)]
        else:
            count = 0
        mean += count
        std += count*count
    total = mean
    mean /= m
    std = np.sqrt(std/m - mean*mean)
    # normalize
    # mean /= N-M-1
    # std /= np.sqrt(N-M-1)
    return(mean, std)

M = 20     # test moduli up tp M
N = 10000  # last block visited

hres1 = {}     # residues modulo m: counts based on model
hres2 = {}     # residues modulo m: counts based on random numbers
fstring = ""   # digit sequence based on model
rstring = ""   # digit sequence based on random numbers

for k in range(2, N):

    # fastint = (3**k) % (2**k)
    # fastint = ((3**k + k)) % (2**k)
    fastint = ((5**k) >> k) % (2**k)
    fstring = fstring + int_to_binstring(fastint, k)
    randint = random.getrandbits(k)
    rstring = rstring + int_to_binstring(randint, k)

    for m in range(2,M):
      
        res1 = fastint % m
        res2 = randint % m
        key1 = (m, res1)
        key2 = (m, res2)
        if k > m:
            hres1 = hash_update(key1, hres1)
            hres2 = hash_update(key2, hres2)

for m in range(2, M):
    for res in range(0,m):
        key = (m, res)
        if key not in hres1:
            hres1[key] = 0
            print(">>> Missing:", key) 
        if key not in hres2:
            hres2[key] = 0
        print(key, hres1[key], hres2[key])

print("\nresidue classes: summary")
for m in range(2, M):
    (fmean, fstd) = summary(hres1, m, N)
    (rmean, rstd) = summary(hres2, m, N)
    # fmean = rmean, thus nor showing rmean
    print("%4d %8.2f %8.2f %9.2f" % (m, fmean, fstd, rstd))

print("\nSubstring counts in digit sequence")
substrings = ('0','1','00','01','10','11',
              '000','001','010','011',
              '100','101','110','111',)
for str in substrings:
    print("%4s %10d %10d" % (str, fstring.count(str), rstring.count(str)))

print(fstring[0:1000])
