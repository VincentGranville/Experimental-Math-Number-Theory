import random  

n = 100000
L = 2*n
H = int(0.04*n)  

sprod = '1'
random.seed(334)

for k in range(n+1):
    u = random.uniform(0, 1)
    if u < 0.00004:
        sprod += '1'
    else: 
        sprod += '0'

scnt1 = sprod.count('1')  # number of 1 in seed
prod = int(sprod, 2)      # seed

arr_count1 = []
arr_count0 = []
arr_count1.append(scnt1)
arr_count0.append(n+1-scnt1)

#--- 1. Main

for k in range(1, H+1): 

    prod = prod*prod 
    pstri = bin(prod)
    stri = pstri[0: L+2]  
    prod = int(stri, 2)
    lsr = len(stri)
    stri = stri[2: len(stri)]
    estri = stri[0:n]   # leftmost n digits
    ecnt0 = estri.count('0')
    ecnt1 = estri.count('1')
    arr_count1.append(ecnt1)
    arr_count0.append(ecnt0)
    print("%3d %3d %3d" %(k, ecnt0, ecnt1)) 

#--- 2. Create the plots

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8

xvalues = np.arange(1, H+2, 1)
arr_color = []

for k in range(H+1):
    if k%2 == 0:
        arr_color.append((1, 0, 0)) # red
    else:
        arr_color.append((0, 0, 1)) # blue 

plt.scatter(xvalues[0:H],arr_count1[0:H],s=0.01,c=arr_color[0:H]) 
plt.axhline(y=n/2,color='red',linestyle='--',linewidth=0.4,dashes=(5,10))
plt.xlim([0, H])
plt.ylim([0.3*n, 0.505*n])
plt.show()

print("\nNumber of 1 in seed:", scnt1)
