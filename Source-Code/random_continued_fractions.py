import numpy as np
np.random.seed(699175)

def fibo(N, step, mode):  ## a, b):
    arr_y = [] 
    arr_y_old = []
    xvals = []

    y = 1  # y = y0 = x1/x0 
    sum = np.log(y)
    sum_y = y
    sum_v = 0
    old_y = 0

    for k in range(1, N):

        if mode == 'uniform':
            v = np.random.uniform(0, 1)
        elif mode == 'bernoulli': 
            v = np.random.randint(0,2)
        elif mode == 'periodic':
            v = 1 + (k % 5)
        sum_v += v

        y = abs(1 + v/y)
        sum_y += y
        sum += np.log(y)
        if k % step == 0 and k > N/2:
            print("%7d %9.7f %9.7f %9.7f" % (k, sum/k, sum_y/k, sum_v/k))
            arr_y.append(y)
            arr_y_old.append(old_y)
            xvals.append(k)
        old_y = y
    return(arr_y, arr_y_old, xvals)


#--- main

N = 600000
step = 1
mode = 'bernoulli' # options: 'uniform', 'periodic', 'bernoulli'
arr_y, arr_y_old, xvals = fibo(N, step, mode)


#--- plot results

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
#plt.rcParams['axes.facecolor'] = 'black'

plt.scatter(xvals, arr_y, s = 0.2)
plt.show()
plt.scatter(arr_y, arr_y_old, s=0.2)
plt.show()


#--- print summary stats

if mode in ('bernoulli', 'periodic'):
    hash_y = {}
    for k in range(len(arr_y)):
        y = arr_y[k]
        y_old = arr_y_old[k]
        key = (y, y_old)
        if key in hash_y:
            hash_y[key]+= 1
        else:
            hash_y[key] = 1
    for key in hash_y: 
        print(key, hash_y[key])
