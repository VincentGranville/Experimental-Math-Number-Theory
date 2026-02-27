import numpy as np
import gmpy2

ctx = gmpy2.get_context()  
ctx.precision = 10000 
ndigits = ctx.precision  
nmatch = 0
newdigits = ''
hash_newdigits = {}
hash_pairs = {}

p = gmpy2.mpz(1)
q = gmpy2.mpz(2)   
mu = 4 
nu = 3 
N = ndigits // nu 

lim = 2**nu - 1 - gmpy2.sqrt(4**nu - 2**(nu+1))  
str_lim = gmpy2.digits(lim, 2)[0]
str_lim = str_lim[0:ndigits]

def update_hash(hash, key, count):
    if key in hash:
        hash[key] += count
    else:
        hash[key] = count
    return()

def count_matching_prefix_chars(str1, str2):
    count = 0
    for char1, char2 in zip(str1, str2):
        if char1 == char2:
            count += 1
        else:
            break
    return(count)

def phi(p, q, nu):
    while p * 2**nu > q:  
         p = p // 2
    return(p)

#--- 1. Main
 
for k in range(N):

    p = (p+q)**2  
    q = 2**mu *q**2 
   
    str_p = gmpy2.digits(p, 2)
    str_p = str_p[0:ndigits]
    p = gmpy2.mpz(str_p, 2)

    str_q = gmpy2.digits(q, 2)
    str_q = str_q[0:ndigits]
    q = gmpy2.mpz(str_q, 2)

    old_p = p
    p = phi(p, q, nu)

    old_nmatch = nmatch
    nmatch = count_matching_prefix_chars(str_p, str_lim)
    old_newdigits = newdigits
    newdigits = str_p[old_nmatch:nmatch]  
    pair = (old_newdigits, newdigits)
    update_hash(hash_pairs, pair, 1)
    x = p/q
    update_hash(hash_newdigits, newdigits, 1)
    if k % 1000 == 0:
        print("New digits:", k, nmatch, newdigits) 

#--- 2. Summary results

print("\n\n")
hash_newdigits = dict(sorted(hash_newdigits.items(), key=lambda item: item[1], reverse=True))
c1 = 0   # counts digits equal to 1
c0 = 0   # counts digits equal to 0
for key in hash_newdigits:
    noccur = hash_newdigits[key]
    c1 += noccur * key.count('1')
    c0 += noccur * key.count('0')
    print("New digits hash summary:",c0, c1, noccur, key)

print("\n\n")
hash_pairs = dict(sorted(hash_pairs.items(), key=lambda item: item[1], reverse=True))
for pair in hash_pairs:
    cnt = hash_pairs[pair]
    if cnt > 5:
        print("Pair:", cnt, pair)
