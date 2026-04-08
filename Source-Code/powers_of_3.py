import numpy as np
import gmpy2
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 'x-small'

ndigits = 50 # maximum possible match
ctx = gmpy2.get_context()  
ctx.precision = ndigits  
N = 200000  
string_tests = ('0', '1', '00', '01', '10', '11', '0000', '1111')

power3 = gmpy2.mpz(1)
isqrt2 = gmpy2.isqrt(2 * 2**(2*ndigits))
rescaled2 = isqrt2/ 2**int(gmpy2.log2(isqrt2))
sqrt2 = bin(isqrt2)[2:]

def update_hash(hash, key, count):
    if key in hash:
        hash[key] += count
    else:
        hash[key] = count 
    return(hash)

def match_strings(bin3, sqrt2):
    # match = number of consecutive identical chars on the left
    match = 0
    Flag = True
    while match < min(len(bin3)-1, len(sqrt2)-1) and Flag:
        if bin3[match] == sqrt2[match]:
            match += 1
        else:
            Flag = False
    return(match)

arr_match = []
arr_nu = []
hash_strings = {}
hash_match = {}
q = 3
lg23 = np.log2(q)
lg2 = np.log(2)
lg3 = np.log(q)
e2 = 2**ndigits


for k in range(N):

    rescaled3 = gmpy2.exp(k*lg3 - int(k*lg23)*lg2)
    arr_nu.append(float(rescaled3))
    power3 = gmpy2.mpz(e2 * rescaled3)
    bin3 = bin(power3)[2:]
    for string in string_tests:
        if bin3[1:len(string)+1] == string:
            update_hash(hash_strings, string, 1)
    match = match_strings(bin3, sqrt2)  
    if match not in hash_match:
        hash_match[match] = k
    arr_match.append(match)

#--- Show frequency of some leading digit strings

print()
for string in hash_strings:
    count = hash_strings[string]
    print("String frequency %5s: %8.6f" % (string, count/N))  

#--- Show phi(kappa)

print()
for kappa in hash_match:
    print("phi(%3d) = %6d" %(kappa, hash_match[kappa]))

#--- Plot convergence of median to sqrt(2)

arr_median_k = []
arr_median_idx = []
arr_median_match = []
arr_median_value = []
arr_median_color = []

for k in range(len(arr_nu)//20, len(arr_nu), 100):
    # compute nu((k-1)/2)
    arr = np.array(arr_nu[0: k])
    median_idx = np.argpartition(arr, len(arr) // 2)[len(arr) // 2]
    match = arr_match[median_idx]
    if match == 18:
        arr_median_color.append('lime')
    elif match == 17:
        arr_median_color.append('yellow')
    elif match == 16:
        arr_median_color.append('magenta')
    else:
        arr_median_color.append('cyan')
    arr_median_k.append(k)
    arr_median_idx.append(median_idx)
    arr_median_match.append(match)
    arr_median_value.append(arr[median_idx])
    if k % 10000 == 1:
        print("Convergence: %6d %6d %8.5f %8.5f" 
                %(k, median_idx, arr[median_idx], arr_match[median_idx]))
 
#--- Plots

plt.rcParams['axes.facecolor'] = 'black'
plt.scatter(arr_median_k, arr_median_value,s=0.6, c='yellow', linewidths=0) 
plt.axhline(y=np.sqrt(2), color='r', linewidth = 0.6)
plt.show()
plt.scatter(arr_median_k, arr_median_idx, s=0.6, linewidths=0, c = arr_median_color)
plt.grid(color='gray')
plt.show()



