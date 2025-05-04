# nt_digit_transform.py: y2 = y(1-y) 

import numpy as np
import random
import gmpy2

random.seed(410)

def randomize(n, q):

    # q is the number of digits set to 1, among the n digits

    u = random.sample(range(1,n-1), q-1)
    v = np.zeros(n-1)

    for k in range(q-1):
        index = u[k]
        v[index] = 1

    stri = [str(int(bit)) for bit in v]
    stri = "".join(stri)
    stri = '1' + stri
    y = gmpy2.mpz(stri, 2)

    return(y, stri)

#--- 1. Main

n = 300 
ctx = gmpy2.get_context()  
ctx.precision = 4*n 
samples = 2000 

arr_dc = np.zeros(n+1)
arr_min = np.zeros(n)
arr_x = []
arr_y = []

low = 1 
high = n-1 
for q in range(low, high):
    if q % 10 == 0:
        print("q=",q)
    min_dc = n+1

    for k in range(samples):
        (y, stri) = randomize(n, q)
        y2 = (2**n * y*(2**(n) - y))
        stri2 = bin(y2)[2:n+2]
        dc = stri2.count('1')
        if dc < min_dc:
            min_dc=dc
            min_stri2 = stri2
            min_stri = stri
        arr_dc[dc]+=1
        arr_x.append(q)
        arr_y.append(dc)

    arr_min[q] = min_dc


#--- 2. Main plot

arr_x = np.array(arr_x)/n
arr_y = np.array(arr_y)/n

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.facecolor'] = 'black'

plt.scatter(arr_x, arr_y, s=0.2, c='bisque', alpha = 0.6)
plt.plot((0, 1), (0, 1), c='red', linewidth=0.6) ##, marker = 'o')
plt.plot((0, 1), (0.5, 0.5), c='red', linewidth=0.6) ##, marker = 'o')
plt.ylim((0.00, 1.00))
plt.xlim((0.00, 1.00))
plt.show()


#--- 3. Other plots

arr_dc_cdf = np.zeros(n)
arr_dc_cdf[0] = arr_dc[0]
for k in range(1, n):
    arr_dc_cdf[k] = arr_dc[k] + arr_dc_cdf[k-1]

xval = range(low, high)
plt.plot(xval, arr_min[low:high])
plt.show()
plt.plot(xval, arr_dc_cdf[low:high])
plt.show()
