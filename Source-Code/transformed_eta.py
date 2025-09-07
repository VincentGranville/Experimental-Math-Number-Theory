from scipy.optimize import fsolve
import numpy as np

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8

#--- 1. Precomputations

def pre_compute(swap_primes, del_primes, N, mode):

    p_alog = np.zeros(N)
    p_arr_omega = np.zeros(N)
    p_arr_delta = np.zeros(N)
    flag = 1 

    for k in range(1, N):

        if k % 100 == 0:
            print("pre_compute:", k)
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

        p_alog[k] = np.log(new_k)
        p_arr_omega[k] = new_k
        p_arr_delta[k] = flag
        flag = -flag

        for p in del_primes:
            if k % p == 0:
                p_arr_delta[k] = 0

    if mode == 'sorted': 
        ranks = np.abs(p_arr_omega).argsort().argsort()
        alog = np.zeros(N)
        arr_omega = np.zeros(N)
        arr_delta = np.zeros(N)
        for k in range(len(arr_omega)):
            rank = ranks[k]
            alog[rank] = p_alog[k]
            arr_omega[rank] = p_arr_omega[k]
            arr_delta[rank] = p_arr_delta[k]
        return(alog, arr_omega, arr_delta)
    else: 
        return(p_alog, p_arr_omega, p_arr_delta)


def eta2(s, params):  

    # t is imaginary part
    sigma = s[0]
    t = s[1]
    n = params[0]
    alog = params[1]
    arr_omega  = params[2]
    arr_delta = params[3]
    vexp = np.exp(t*alog[1:n]*1j) 
    vpow = arr_delta[1:n] * arr_omega[1:n]**(-sigma) 
    sum = np.dot(vexp, vpow)
    norm = np.linalg.norm(sum)
    return(norm)

#--- 2. Removing primes in eta function

sigma = 0.5 
t = 25.010858 # 0.5 + it is 3rd root of zeta 
s = [sigma, t]
n0 = 4
n1 = 3500
N = 500001
swap_primes = {} 
mode = 'unsorted' # options: 'sorted', 'unsorted' 

del_primes1 = (3, 5, 7) 
m1 = len(del_primes1) + 1
alog, arr_omega, arr_delta = pre_compute(swap_primes, del_primes1, N, mode)
params1 = [0, alog, arr_omega, arr_delta]

del_primes2 = () 
m2 = len(del_primes2) + 1
alog, arr_omega, arr_delta = pre_compute(swap_primes, del_primes2, N, mode)
params2 = [0, alog, arr_omega, arr_delta]

arr_x = []
arr_y1 = []
arr_y2 = []
for n in range(n0, n1):
    params1[0] = n
    params2[0] = n
    arr_x.append(n-1)
    arr_y1.append(eta2(s, params1)) 
    arr_y2.append(eta2(s, params2)) 

plt.scatter(arr_x, arr_y1, c='green', s = 0.1, label="m = %d" %m1)  
plt.plot(arr_x, arr_y2, c='red', linewidth = 0.6, label="m = %d" %m2)
plt.axhline(y=0.0, color='black', linestyle='dotted', linewidth = 0.3) 
plt.legend(loc="upper right")
plt.show()

#--- 3. Changing primes in eta function

def eta(s, params): 
    # s[0] = real part, s[1] = imaginary part
    n = params[0]
    alog = params[1]
    arr_omega  = params[2]
    arr_delta = params[3]
    vexp = np.exp(s[1]*alog[1:n]*1j) 
    vpow = arr_delta[1:n] * arr_omega[1:n]**(-s[0]) 
    sum = np.dot(vexp, vpow)
    return [sum.real, sum.imag]

swap_primes = { 
                  3: 4.051, 
                  5: 7.916, 
                  7: 9.114 
              } 
del_primes = () 
mode = 'unsorted' # options: 'sorted', 'unsorted' 
alog, arr_omega, arr_delta = pre_compute(swap_primes, del_primes, N, mode)
arr_x = []
arr_y = []
arr_z = []

sigma = 0.5
t1 = 14.13472514
t3 = 25.01085758
s1 = [sigma, t1]
s3 = [sigma, t3]
params = [0, alog, arr_omega, arr_delta]
for n in range(400, 300000, 1000):
    arr_x.append(n-1)
    params[0] = n
    arr_y.append(eta2(s1, params)) 
    arr_z.append(eta2(s3, params))

plt.plot(arr_x, arr_y, c='green', linewidth = 0.8, label="s = s1")  
plt.plot(arr_x, arr_z, c='blue', linewidth = 0.8, label="s = s3")
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

params = [n, alog, arr_omega, arr_delta]
for real in np.arange(0.45, 1.05, 0.001): 
    print("Real: %7.4f" %(real))
    for imag in np.arange(10.00, 40.01, 0.1):
        init = np.array([real, imag])
        root = fsolve(eta, init, params)
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

plt.scatter(aimag, areal, s = 0.2, c=acolor)
plt.show()
print()
print("Found: left: %5.3f right: %5.3f" % (found1/nfc, found2/nfc))

#--- 5. Replace all primes by random numbers

from sympy import primerange, isprime

def build_swap_primes(nprimes, noise):

    arr_primes = list(primerange(3, nprimes))
    swap_primes = {}
    seed = 503
    np.random.seed(seed)
    N_max = len(arr_primes)
    for k in range(N_max):
        prime = arr_primes[k]
        u = np.random.uniform(-noise, noise)
        swap_primes[prime] = prime + u*(k**0.25) ##### k**0.25
    print("N_max: ", N_max)
    return(swap_primes, N_max)

del_primes = ()
nprimes = 300000
sigma = 0.5 
t = 25.010858 # 0.5 + it is 3rd root of zeta
s = [sigma, t]

noise1 = 0.005 
swap_primes1, N_max = build_swap_primes(nprimes, noise1)
mode = 'sorted' # options: 'sorted', 'unsorted'
N = N_max  
n0 = 4
n1 = N_max

alog, arr_omega, arr_delta = pre_compute(swap_primes1, del_primes, N, mode)
params1 = [0, alog, arr_omega, arr_delta]

noise2 = 0.000 
swap_primes2, N_max = build_swap_primes(nprimes, noise2)
alog, arr_omega, arr_delta = pre_compute(swap_primes2, del_primes, N, mode)
params2 = [0, alog, arr_omega, arr_delta]

arr_x = []
arr_y1 = []
arr_y2 = []
for n in range(n0, n1):
    params1[0] = n
    params2[0] = n
    arr_x.append(n)
    arr_y1.append(eta2(s, params1))
    arr_y2.append(eta2(s, params2))

plt.scatter(arr_x, arr_y1, c='green', label="noise = %5.3f" %noise1, s=0.1) ### linewidth = 0.3)  
plt.plot(arr_x, arr_y2, c='red', label="noise = %5.3f" %noise2, linewidth = 0.8)  
plt.axhline(y=0.0, color='black', linestyle='dotted', linewidth = 0.3) 
plt.legend(loc="upper right")
plt.show()
