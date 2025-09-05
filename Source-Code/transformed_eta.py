from scipy.optimize import fsolve, minimize, Bounds
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8

#--- 1. Precomputations

def pre_compute(swap_primes, del_primes):

    alog = np.zeros(N)
    arr_omega = np.zeros(N)
    arr_delta = np.zeros(N)
    flag = 1 

    for k in range(1, N):

        hash_pval = {}
        new_k = k
        for p in swap_primes:
            pval = 0
            while new_k % p == 0:
                new_k = new_k // p
                pval += 1
            hash_pval[p] = pval

        for p in hash_pval:
            new_k *= (swap_primes[p])**hash_pval[p]

        alog[k] = np.log(new_k)
        arr_omega[k] = new_k
        arr_delta[k] = flag
        flag = -flag

        for p in del_primes:
            if k % p == 0:
                arr_delta[k] = 0

    return(alog, arr_omega, arr_delta)

def eta2(sigma, t):

    # t is imaginary part
    vexp = np.exp(t*alog[1:n]*1j) 
    vpow = arr_delta[1:n] * arr_omega[1:n]**(-sigma) 
    sum = np.dot(vexp, vpow)
    norm = np.linalg.norm(sum)
    return(norm)

#--- 2. Removing primes in eta function

sigma = 0.5 
t = 25.010858 # 0.5 + it is 3rd root of zeta 
n0 = 4
n1 = 3500
N = 500001
swap_primes = {} 

del_primes = (3, 5, 7) 
m = len(del_primes) + 1
alog, arr_omega, arr_delta = pre_compute(swap_primes, del_primes)
for k in range(1, 40):
    if arr_delta[k] != 0:
        print("%3d %12.7f" %(k, arr_omega[k]*arr_delta[k]))

arr_x = []
arr_y = []
for n in range(n0, n1):
    arr_x.append(n-1)
    arr_y.append(eta2(sigma, t))

del_primes = ()
alog, arr_omega, arr_delta = pre_compute(swap_primes, del_primes)
arr_z = []
for n in range(n0, n1):
    arr_z.append(eta2(sigma, t))

plt.plot(arr_x, arr_y, c='green', linewidth = 0.3, label="m = %d" %m)  
plt.plot(arr_x, arr_z, c='red', linewidth = 0.6, label="m = 1")
plt.axhline(y=0.0, color='black', linestyle='dotted', linewidth = 0.3) 
plt.legend(loc="upper right")
plt.show()

#--- 3. Changing primes in eta function

def eta(s, n):
    # s[0] = real part, s[1] = imaginary part
    vexp = np.exp(s[1]*alog[1:n]*1j) 
    vpow = arr_delta[1:n] * arr_omega[1:n]**(-s[0]) 
    sum = np.dot(vexp, vpow)
    return [sum.real, sum.imag]

swap_primes = { 
                  3: 4.051, # np.pi, 
                  5: 7.916, # np.log(131), 
                  7: 9.114 # np.exp(2)
              } 
del_primes = () 
alog, arr_omega, arr_delta = pre_compute(swap_primes, del_primes)
arr_x = []
arr_y = []
arr_z = []

sigma = 0.5
t1 = 14.13472514
t3 = 25.01085758
for n in range(400, 100000, 1000):
    arr_x.append(n-1)
    arr_y.append(eta2(sigma, t1))
    arr_z.append(eta2(sigma, t3))

plt.plot(arr_x, arr_y, c='green', linewidth = 0.3, label="s = s1")  
plt.plot(arr_x, arr_z, c='red', linewidth = 0.6, label="s = s3")
plt.axhline(y=0.0, color='black', linestyle='dotted', linewidth = 0.3) 
plt.legend(loc="upper right")
plt.show()

#--- 4. Basins of attraction

areal = []
aimag = []
acolor = []
n = 499
found1 = 0
found2 = 0
nfc = 0
eps = 0.05

for real in np.arange(0.45, 1.05, 0.001):
    print("Real: %7.4f" %(real))
    for imag in np.arange(10.00, 40.01, 0.1):
        init = np.array([real, imag])
        root = fsolve(eta, init, n)
        areal.append(real)
        aimag.append(imag)
        sigma = min(1, root[0])
        nfc += 1
        if 0 < abs(root[0] - 0.5) < eps:
            found1 += 1
            acolor.append('red')
        elif 0 < abs(root[0] - 1.0) < eps:
            found2 += 1
            acolor.append('green')
        elif 0 < root[0] < 1:
            acolor.append('blue')
        else:
            acolor.append('gold')

plt.scatter(areal, aimag, s = 0.2, c=acolor)
plt.show()
print()
print("Found: left: %5.3f right: %5.3f" % (found1/nfc, found2/nfc))


