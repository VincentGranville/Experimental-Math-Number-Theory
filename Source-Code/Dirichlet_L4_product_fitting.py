# DirichletL4_EulerProduct.py
# On WolframAlpha: DirichletL[4,2,s], s = sigma + it
#     gives the Dirichlet L-function  for the Dirichlet character  with modulus k and index j.
#
# References:
#     https://www.maths.nottingham.ac.uk/plp/pmzcw/download/fnt_chap4.pdf
#     https://mpmath.org/doc/current/functions/zeta.html
#     f(s) = dirichlet(s, [0, 1, 0, -1]) in MPmath
#     https://mathoverflow.net/questions/439018/zeros-of-dirichlet-function-ls-chi-4/439019#439019

import matplotlib.pyplot as plt
import matplotlib as mpl
import mpmath
import numpy as np
from primePy import primes
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import warnings
warnings.filterwarnings("ignore")

#--- [1] create tables of prime numbers

m = 1000000  # primes up to m included in Euler product
aprimes  = []  

for k in range(m):
    if k % 100000 == 0:
        print("Creating prime table up to p <=", k)
    if primes.check(k) and k > 2:
        aprimes.append(k)


#--- [2] Euler product

#--- [2.1] Main function

def L4_Euler_prod(mode = 'L4', sigma = 1.00, t = 0.00):

    L4 = mpmath.dirichlet(complex(sigma,t), [0, 1, 0, -1]) 
    print("\nMPmath lib.: L4(%8.5f + %8.5f i) = %8.5f + %8.5f i" 
          % (sigma, t,L4.real,L4.imag))

    prod = 1.0
    sum_chi4 = 0 
    sum_delta = 0
    run_chi4 = 0
    old_chi4 = 0
    DLseries = 0
    flag = 1

    aprod = [] 
    adelta = []
    asum_delta = []
    achi4 = []
    arun_chi4 = []
    asum_chi4 = []

    x1 = []
    x2 = []
    error1 = []
    error2 = []
    seed = 116  # try 103, 105, 116 & start = 2000 (for mode = 'rn')
    np.random.seed(seed)
    eps = 0.000000001

    for k in range(len(aprimes)): 
        
        if mode == 'L4':
            condition  = (aprimes[k] % 4 == 1)  
        elif mode == 'Q2':
            condition  = (k % 2 == 0)           
        elif mode == 'rn':
            condition = (np.random.uniform(0,1) < 0.5)
        
        if condition:   
            chi4 = 1
        else:
            chi4 = -1

        sum_chi4 += chi4 
        achi4.append(chi4)
        omega = 1.00 # try 1.00, sigma or 1.10
        # if omega > 1, asum_chi4[n] --> 0 as n --> infty
        # asum_chi4.append(sum_chi4/aprimes[k]**omega)
        asum_chi4.append(sum_chi4/(k+1)**omega) 
        # asum_chi4.append(sum_chi4/(k+1)*(np.log(k+2))) 

        if chi4 == old_chi4:
            run_chi4 += chi4
        else:
            run_chi4 = chi4
        old_chi4 = chi4
        arun_chi4.append(run_chi4) 
 
        factor = 1 - chi4 * mpmath.power(aprimes[k], -complex(sigma,t))
        prod *= factor
        aprod.append(1/prod)    

        term = mpmath.power(2*k+1, -complex(sigma,t))
        DLseries += flag*term
        flag = -flag

    limit = -eps + 1/prod   # full finite product (approx. of the limit)
    if mode == 'L4':
        limit = L4   # use exact value instead (infinite product if it converges)
  
    for k in range(len(aprimes)): 

        delta = (aprod[k] - limit).real  # use real part
        adelta.append(delta)
        sum_delta += delta
        asum_delta.append(sum_delta)
        chi4 = achi4[k]
       
        if chi4 == 1: 
            x1.append(k)
            error1.append(delta)
        elif chi4== -1: 
            x2.append(k)
            error2.append(delta)

    print("Dirichlet L: DL(%8.5f + %8.5f i) = %8.5f + %8.5f i" 
        % (sigma, t, DLseries.real, DLseries.imag))
    print("Euler Prod.: %s(%8.5f + %8.5f i) = %8.5f + %8.5f i\n" 
        % (mode, sigma, t, limit.real, limit.imag))

    adelta = np.array(adelta)
    aprod = np.array(aprod)
    asum_chi4 = np.array(asum_chi4)  
    asum_delta = np.array(asum_delta) 
    error1 = np.array(error1)
    error2 = np.array(error2)

    return(limit.real, x1, x2, error1, error2, aprod, adelta, asum_delta, 
             arun_chi4, asum_chi4)

#--- [2.2] Main part
    
mode = 'L4'  # options: 'L4', 'Q2', 'rn' (random chi4)
(prod, x1, x2, error1, error2, aprod, adelta, asum_delta, arun_chi4, 
     asum_chi4) = L4_Euler_prod(mode, sigma = 0.90, t = 0.00)


#--- [3] Plots (delta is Euler product, minus its limit)

mpl.rcParams['axes.linewidth'] = 0.3
plt.rcParams['xtick.labelsize'] = 7
plt.rcParams['ytick.labelsize'] = 7

#- [3.1] Plot delta and cumulated chi4

x = np.arange(0, len(aprod), 1)

# offset < len(aprimes), used to enhance visualizations 
offset = int(0.02 * len(aprimes))  

# y1 = aprod / prod
# plt.plot(x[offset:], y1[offset:], linewidth = 0.1)
# plt.show()

y2 = adelta 
plt.subplot(2,2,1)
plt.plot(x[offset:], y2[offset:], marker=',', markersize=0.1, 
   linestyle='None', c='red')

y3 = asum_chi4 
plt.subplot(2,2,2)
plt.plot(x[offset:], y3[offset:], marker=',', markersize=0.1, 
   linestyle='None', c='red')

#- [3.2] Denoising L4, curve fitting

def objective(x, a, b, c):

    # try c = 0 (actual limit)
    value = c + a/np.sqrt(x) +  b/np.sqrt(x*np.log(x)) 
    return value

def model_fit(x, y2, y3, start, offset, n_max):
   
    for k in range(n_max):

        n = int(len(y2) * (k + 1) / n_max) - start  
        stdn_y2 = np.std(y2[start:n])
        stdn_y3 = np.std(y3[start:n])
        r = stdn_y3 / stdn_y2
        
        # note: y3 / r ~ mu / sqrt(x)  [chaotic part]
        mu = y3[n] * np.sqrt(n)  # tend to a constant ?
        y4 = y2 * r - y3
        y4_fit = []
        err = -1

        if min(y4[start:]) > 0: 
            popt, pcov = curve_fit(objective, x[start:n], y4[start:n], 
                         p0=[1, 1, 0], maxfev=5000)
            [a, b, c] = popt
            y4_fit = objective(x, a, b, c)
            err = r2_score(y4[offset:], y4_fit[offset:])
            print("n = %7d  mu =%6.2f  c =%6.2f  a =%5.2f  b =%5.2f  r =%6.3f  err =%6.3f"
                        %(n, mu, c, a, b, r, err))

    return(y4, y4_fit, err, n)

n_max = 10   # testing n_max values of n, equally spaced
start = 20   # use Euler products with at least 'start' factors 
if mode == 'rn':
    start = 1000
if start > 0.5 * offset:
    print("Warning: 'start' reduced to 0.5 * offset")
    start = int(0.5 * offset)
(y4, y4_fit, err, n) = model_fit(x, y2, y3, start, offset, n_max)
ns = np.sqrt(n)

if err != -1:
    plt.subplot(2,2,3)
    plt.plot(x[offset:], ns*y4[offset:], marker=',', markersize=0.1, 
       linestyle='None', c='orange')
    plt.plot(x[offset:], ns*y4_fit[offset:], linewidth = 0.2, c='black')
else:
    print("Can't fit: some y4 <= 0 (try different seed or increase 'start')")

#--- [3.3] Plot integral of delta

y5 = asum_delta
plt.subplot(2,2,4)
plt.plot(x[offset:], y5[offset:], linewidth = 0.4, c='red')
plt.show()


#--- [4] Quantum derivative

#- [4.1] Function to differentiated: delta, here broken down into 2 legs

plt.subplot(1,2,1)
shift = 0.001
plt.plot(x1[offset:], error1[offset:], marker=',', markersize=0.1, 
   linestyle='None', alpha = 1.0, c='red')
plt.plot(x2[offset:], shift + error2[offset:], marker=',', markersize=0.1, 
   linestyle='None', alpha = 0.2, c='orange')

#- [4.2] Quantum derivative

def d_error(arr_error):

    diff_error = []  # discrete derivative of the error
    positives = 0
    negatives = 0
    for k in range(len(arr_error)-1):
        diff_error.append(arr_error[k+1] - arr_error[k])
        if arr_error[k+1] - arr_error[k] > 0: 
            positives +=1
        else:
            negatives += 1
    return(diff_error, positives, negatives)

(diff_error1, positives1, negatives1) = d_error(error1)
(diff_error2, positives2, negatives2) = d_error(error2)
ymin = 0.5 * float(min(min(diff_error1[offset:]), min(diff_error1[offset:])))
ymax = 0.5 * float(max(max(diff_error1[offset:]), max(diff_error2[offset:])))

plt.subplot(1,2,2)
plt.ylim(ymin, ymax)
plt.plot(x1[offset:len(x1)-1], diff_error1[offset:len(x1)-1], marker=',', markersize=0.1, 
         linestyle='None', alpha=0.8, c = 'red')
plt.plot(x2[offset:len(x2)-1], diff_error2[offset:len(x2)-1], marker=',', markersize=0.1, 
         linestyle='None', alpha=0.8, c = 'orange')
plt.show()

print("\nError 1: positives1: %8d  negatives1: %8d" % (positives1, negatives1)) 
print("Error 2: positives2: %8d  negatives2: %8d" % (positives2, negatives2)) 

