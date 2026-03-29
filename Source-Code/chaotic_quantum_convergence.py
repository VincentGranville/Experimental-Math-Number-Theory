import numpy as np
import gmpy2
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 'x-small'

kmax = 5000
ctx = gmpy2.get_context()  
ctx.precision = kmax 
ndigits = ctx.precision  
z = gmpy2.log(2)

N = 5000
mode = 'quantum'  # options: 'quantum' or 'chaotic'

x0 = gmpy2.mpfr(0)
x1 = gmpy2.mpfr(1)
x2 = gmpy2.mpfr(1)
arr_k = []
arr_x = []
arr_col = []

def set_color(k, mode):

    if mode == 'quantum':
        if k % 7 == 0:
            color='red'
        elif k % 7 == 1:
            color='blue'
        elif k % 7 == 2:
            color='green'
        elif k % 7 == 3:
            color='gold'
        elif k % 7 == 4:
            color='orange'
        elif k % 7 == 5:
            color='magenta'
        elif k % 7 == 6:
            color='gray'

    elif mode == 'chaotic':
        if k % 4 == 0:
            color = 'red'
        elif k % 4 == 1:
            color = 'gold'
        elif k % 4 == 2:
            color = 'green'
        elif k % 4 == 3:
            color = 'magenta'
        else:
            color = 'gray'
    return(color)


for k in range(2,N):
    color = set_color(k, mode)
    if mode == 'quantum':
        x = 2*x2 - 16*x1 + 4*x0 
    elif mode == 'chaotic':
        x = abs(x2 - 3*x1)
    v = gmpy2.log(abs(x))/k
    w = gmpy2.exp(v)
    if k % 100 == 0:
        print("rho_n", k, w)   
    arr_k.append(k)
    arr_x.append(w)  
    arr_col.append(color)
    x0 = x1
    x1 = x2
    x2 = x

plt.scatter(arr_k, arr_x, c=arr_col, s=0.05)
plt.show()
