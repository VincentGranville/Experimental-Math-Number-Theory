import gmpy2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import pyplot

def update_hash(hash, key, count):
    if key in hash:
        hash[key] += count
    else:
        hash[key] = count
    return(hash)

def p_adic(k, base):
    # return the largest integer v such that base**v divides k
    # https://www.geeksforgeeks.org/python-program-to-find-whether-a-no-is-power-of-two/
    if base == 2:
        div = k & (~(k - 1))  
        v = len(bin(div)) - 3  
    else: 
        div = 1 
        v = 0
        while k %  div == 0:
            div *= base
            v += 1
        div = div // base 
        v -= 1 
    return(v, div)
      
def initialize(mode, z=1, m=0):

    p_buffer = []
    q_buffer = []

    if 'Exponential' in mode:
        # p/q tends to exp(x)-1-1/x as k -> infty  
        p_buffer.append(0)
        p_buffer.append(0)
        q_buffer.append(0)
        q_buffer.append(z)

    elif mode == 'Linear':
        # p/q tends to 1/(3-SQRT(3)) as k --> infty
        p_buffer.append(0) 
        p_buffer.append(1) 
        q_buffer.append(0) 
        q_buffer.append(2) 

    elif mode == 'ContinuedFractions':
        # p/q tends to (-1+sqrt(5))/2 as k -> infty
        p_buffer.append(0) 
        p_buffer.append(13) 
        q_buffer.append(0) 
        q_buffer.append(21) 

    elif mode == 'Special':
        p_buffer.append(0) 
        p_buffer.append(0)  
        q_buffer.append(1) 
        q_buffer.append(2)  

    elif mode == 'Special2':
        p_buffer.append(0) 
        p_buffer.append(0)  
        q_buffer.append(1) 
        q_buffer.append(2)  

    elif mode == 'Special3':
        p_buffer.append(0) 
        p_buffer.append(1)  
        q_buffer.append(0) 
        q_buffer.append(2)  

    return(p_buffer, q_buffer)


#--- [1] main loop

# parameters
mode =  'Exponential'  # options: 'Exponential', 'Exponential1', 'Exponential2'
                       # 'Linear', 'ContinuedFractions', 'Special', 'Special2', 'Special3'
random = False 
base = 2     # base must be a prime number
n = 45000    # number of iterations
m = 0        # used in Exponential mode (if m=0, p/q tends to exp(z)-1-1/z)  
z = 1        # used in Exponential mode (integer > 0)
seed = 659   # used if random = True
np.random.seed(seed)

# local variables
###arr_delta = []
digits = ""
new_digits = ""
hash1 = {}
hash2 = {}
hash_nu = {}
bprefix = ""
lmax = 0
k_old = 0
zeros = 0
ones = 0

(p_buffer, q_buffer) = initialize(mode, z, m) 

for k in range(2, n+1, 1):

    # p = p[k], p_buffer[1] = p[k-1], p_buffer[2] = p[k-2]
    # q = q[k], q_buffer[1] = q[k-1], q_buffer[2] = q[k-2]

    if mode == 'Exponential':
        p = z*k*p_buffer[1] + 1 
        q = z*k*q_buffer[1]

    elif mode == 'Exponential1' and m > 0:
        p = z*min(k, m)*p_buffer[1] + 1  
        q = z*min(k, m)*q_buffer[1]   

    elif mode == 'Exponential2' and m > 0:
        p = z*((k-1)%m + 1) * p_buffer[1] + 1 
        q = z*((k-1)%m + 1) * q_buffer[1]   

    elif mode == 'Linear': 
        p = 2*p_buffer[1] + 2*p_buffer[0] + 1  
        q = 2*q_buffer[1] + 2*q_buffer[0] 

    elif mode == 'ContinuedFractions':
        p = q_buffer[1]   
        q = p_buffer[1] + q_buffer[1] 

    elif mode == 'Special':
        if (2*p_buffer[1]+ 1)**2 < 2 * q_buffer[0]**2:
            p = 2*p_buffer[1] + 1
        else:
            p = 2*p_buffer[1]
        q = 2*q_buffer[1]

    elif mode == 'Special2':
        p = 2**k * p_buffer[1] + ((3**k - 1)//2) % (2**k)
        q = 2**k * q_buffer[1]

    elif mode == 'Special3':
        p = 2**k * p_buffer[1] + (3**k)  % (2**k)
        q = 2**k * q_buffer[1]

    p_buffer[0] = p_buffer[1]
    p_buffer[1] = p
    q_buffer[0] = q_buffer[1]
    q_buffer[1] = q
    nu_k, div_k = p_adic(k, base)
    nu_p, div_p = p_adic(p, base)
    nu_q, div_q = p_adic(q, base)
    if p > q:
        print("Warning: p >= q (unauthorized)")
        exit()
    l = max(0, nu_q - nu_p)   # length of prefix

    if l > lmax:  

        # process additional digits found
        lmax = l
        prefix = (base**l) * p // q 
        bprefix_old = bprefix   
        l_old = len(bprefix_old)
        bprefix = gmpy2.mpz(prefix).digits(base)
        while len(bprefix) < l:
            bprefix = '0' + bprefix
        if bprefix[0:l_old] != bprefix_old:
            match = 'Fail'
        else:
            match = 'Success'
        new_digits_old = new_digits

        if random:
            new_int = np.random.randint(base**(l-l_old)) 
        else:
            new_int = prefix % (base**(l - l_old)) 

        new_digits = gmpy2.mpz(new_int).digits(base)
        while len(new_digits) < l - l_old:
            new_digits = '0' + new_digits

        delta = k - k_old
        zeros += new_digits.count('0')
        ones  += new_digits.count('1')
        ### arr_delta.append(zeros-ones) #################
        size = len(new_digits)
        print("===>", k, l, size, delta, nu_k, nu_p, nu_q, match, zeros+ones, 
               zeros, ones, ">>", new_digits) 

        digits += new_digits
        k_old = k
        key = (size, k)
        if new_digits not in hash1:
            hash_nu[key] = new_digits
        update_hash(hash1, new_digits, 1)
        key = (new_digits, new_digits_old)
        update_hash(hash2, key, 1)


#--- [2] Output results

#- [2.1] block counts

patterns = () # list of all potential combos of t binary digits
t = 5 
for k in range(0,2**t,1):
    bint = bin(k)
    bint = bint[2:len(bint)]
    while len(bint)<t:
        bint = "0"+bint
    patterns = (*patterns, bint)

print("\nblock counts\n")
hash1_print = {}

for key in hash1:
    klen = len(key)
    hash1_print[(klen, key)] = hash1[key] 

for keyx in sorted(hash1_print):
    key = keyx[1]
    klen = len(key) 
    count = hash1[key]
    if count > 1: 
        print(klen, key, count)

#- [2.2] first occurrence of block

print()
print("first occurrence of block")
old_size = 0
count = 0
for key in sorted(hash_nu):
    size = key[0]
    if size == old_size:
        count = count+1
    else: 
        print()
        count = 1
        old_size = size
    print(count, hash_nu[key], key)

#- [2.3] conditional block counts

print("\nconditional block counts\n")
hash2_print = {}

for key in hash2:
    klenA = len(key[0])
    klenB = len(key[1])
    hash2_print[(klenA, klenB, key)] = hash2[key] 

for keyx in sorted(hash2_print):
    key = keyx[2] 
    old = key[0]
    new = key[1]
    count = hash2[key]
    if count > 5:
        print("(old, new):",key, hash2[key])

#- [2.4] high level summary

print()

print("Zeros vs ones:", zeros, ones)
print("Proportion of zeros:", zeros/(zeros + ones))

sum1 = 0
ndigits1 = min(80,len(digits))

for k in range(0,ndigits1,1):
    sum1 += int(digits[k])/base**(k+1) 

sum2 = 0
ndigits2 = min(80,len(bprefix))

for k in range(0,ndigits2,1):
    sum2 += int(bprefix[k])/base**(k+1)

print("dCheck:", sum1)
print("bCheck:", sum2)
print("Number:", p/q)

if 'Exponential' in mode:
    print("Target:", np.exp(1/z)-(1+1/z))
elif mode == 'Linear':
    print("Target:", 1/(3-np.sqrt(3)))
elif mode == 'ContinuedFractions':
    print("Target:", (-1+np.sqrt(5))/2)
elif mode == 'Special':
    print("Target:", np.sqrt(2)/4)

print("Digits per n:", (zeros+ones)/n)




