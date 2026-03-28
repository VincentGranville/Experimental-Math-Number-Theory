import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 'x-small'

def countSetBits(n):
    count = 0
    while(n):
        count += n & 1
        n >>= 1
    return(count)

sum = 0
arr_k = []
arr_sum = []

for k in range(0,8000000):
    d = countSetBits(k)
    sum += d
    arr_k.append(k)
    if k < 2: 
        arr_sum.append(0)
    else:
        arr_sum.append(sum/(k*np.log2(k)))
    if k % 10 == 0:
        print("sss", k, d, sum)

plt.plot(arr_k, arr_sum, linewidth = 0.2)
plt.show()
