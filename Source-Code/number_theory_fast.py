n = 10000
L = 2*n
H = int(1.05*n)  

def assign_color(k):

    if k%12 == 0:
        color = 'lime'
    elif k%12 == 2:
        color = 'white' 
    elif k%12 == 4:
        color = 'black'
    elif k%12 == 6:
        color = 'cyan'
    elif k%12 == 8:
        color = 'gold'
    elif k%12 == 10:
        color = 'orange'
    elif k%12 == 1:
        color = 'magenta'
    elif k%12 == 3:
        color = 'paleturquoise'
    elif k%12 == 5: 
        color = 'khaki'
    elif k%12 == 7:
        color = 'yellow'
    elif k%12 == 9:
        color = 'mistyrose'
    elif k%12 == 11:
        color = 'orangered'
    return(color)


#--- 1. Main

import gmpy2
import numpy as np

prod = (2**n + 1)               # seed 
ctx = gmpy2.get_context()  
ctx.precision = 2*n       # precision set to 2n bits
prod = gmpy2.mpz(2**n + 1) 

arr_count1 = []
arr_colors = []
xvalues = []
ecnt1 = -1

for k in range(1, H+1): 

    prod = prod*prod
    pstri = bin(int(prod))
    stri = pstri[0: L+2]  
    prod = int(stri, 2)
    prod = gmpy2.mpz(prod) 

    if k>0:  ## 0.8*n: 
        stri = stri[2:]
        if k == n:
            e_approx = stri
        estri = stri[0:n]   # leftmost n digits
        ecnt1 = estri.count('1')
        arr_count1.append(ecnt1)
        color = assign_color(k)
        arr_colors.append(color)
        xvalues.append(k)
 
    if k%1000 == 0:
        print("%3d %3d" %(k, ecnt1))


#--- 2. Fast computation of binary digits of e

from mpmath import mp

# Set precision for n binary digits
mp.dps = int(n*np.log2(10))
e_value = mp.e  # Get e in decimal

# Convert to binary
e_binary = bin(int(e_value * (2 ** n)))[2:]

k = 0
print()
while e_approx[k] == e_binary[k]:
    k += 1
print("%d correct digits (n = %d)" %(k, n))


#---3. Create the plots

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.facecolor'] = 'black'

#ax = plt.axes()
#ax.set_facecolor('silver')

plt.scatter(xvalues, arr_count1,s=0.01,c=arr_colors) 
plt.axhline(y=n/2,color='red',linestyle='--',linewidth=0.6,dashes=(5,10))   
plt.axhline(y=n/5, color='black', linestyle='--', linewidth = 0.6, dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))

for k in range(2,15):
    plt.axvline(x=k*n/(k+1),c='gray',linestyle='--',linewidth=0.6,dashes=(5, 10))

plt.xlim([0.65*n, 1.04*n])
plt.ylim([0, 0.52*n]) 
plt.show()
