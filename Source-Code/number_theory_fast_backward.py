import gmpy2
from gmpy2 import mpfr
import colorsys 

#--- 1. Create table of contrasted colors

def hsv_to_rgb(h, s, v):
    return tuple(round(i * 255) for i in colorsys.hsv_to_rgb(h, s, v))

def generate_contrasting_colors(ncolors):
    colors = []
    for i in range(ncolors):
        hue = i / ncolors
        col = hsv_to_rgb(hue, 1.0, 1.0)
        color = (col[0]/255, col[1]/255, col[2]/255)
        colors.append(color)
    return colors

#--- 2. Functions relared to primorials

def update_q(q, k, file):

    q = gmpy2.mpz(k*q)
    iq = 2**gmpy2.floor(gmpy2.log2(q)) # use floor, not mpz (mpz = round)
    f = str(gmpy2.mpfr(q/iq))[0:20]
    file.write(str(k) + "\t" + f + "\n") 
    return(q)

def primorial(kappa, mode="primorial"): 

    # mode = "primorial" --> return p = #kappa with correct precision
    # mode = "factorial" --> return p = kappa! with correct precision

    from primePy import primes
    ctx = gmpy2.get_context() 
    old_precision = ctx.precision 
    ctx.precision = 2*kappa
    q = gmpy2.mpz(1)   

    OUT = open("primorials.txt", "w")
    # values of r = q/2^int(log2 q) distributed in [1, 2] like F(r) = log2 r
    # this is called the reciprocal distribution
    for k in range(2, kappa+1):
        if mode == "primorial": 
            if primes.check(k): 
                q = update_q(q, k, OUT)
        elif type == "factorials":
            q = update_q(q, k, OUT)
    OUT.close()

    iq = int(gmpy2.floor(gmpy2.log2(q)))
    print("Primorial precision: %d bits | min needed: %d bits" %(ctx.precision, iq))  
    print() 
    p = q
    ctx.precision = old_precision
    return(p)

#--- 3. Main function: the backwars iterations

def iexp(n, start, iters, ncolors, z, u, v, truncate):  

    arr_count1 = []
    arr_colors = []
    xvalues = []
    ecnt1 = -1

    pow2 = 2**(start)  
    z = gmpy2.exp(gmpy2.log(z)/pow2)  # z = exp[p^(1/2^start)]

    for k in range(n-start, n-start-iters, -1):

        iz = gmpy2.mpz(gmpy2.mpfr(2**(n+5) * z)) ### why n+5 ??

        if k % u == v: 
            if truncate:
                # strip 1 and first n-k digits (zeros) on the left
                stri = bin(iz)[2+n-k:2*n-k+2+1] 
                # stri = bin(iz)[2+n-2*k:2*n-k+2+1]
                ecnt1 = stri.count('1') * n/len(stri)
            else:
                stri = bin(iz)[2:n+2+1] 
                ecnt1 = stri.count('1')

            arr_count1.append(ecnt1)
            color = colorTable[k % ncolors]
            arr_colors.append(color)
            xvalues.append(k)

        if k%1000 == 0:
            print(k, ecnt1)            
        z = gmpy2.sqrt(z)

    return(arr_count1, arr_colors, xvalues)

#--- 4. Function to create reverse seed z

# initialize seed z

def initialize_reverse_seed(seed_type):

    if seed_type == "primorial":
        kappa = 30 # try 3, 10, 15, 30, 300, 3000 
        p = primorial(kappa)
        iplog = int(gmpy2.log2(p))
        p = gmpy2.mpfr(p/(2**iplog))
        # try replacing p by -p
        z = gmpy2.exp(p)

    elif seed_type == "random":
        import numpy as np
        np_seed = 6696
        stri =""
        np.random.seed(np_seed)
        for k in range(2*n+1):
            d = np.random.randint(2)
            stri += str(d)
        p = gmpy2.mpz(int(stri, 2)) ###
        iplog = int(gmpy2.log2(p))
        p = gmpy2.mpfr(p/(2**iplog))
        z = gmpy2.exp(p)

    elif seed_type == "integer":
        # also try -1 (backward/forward algo show different paths)
        z = gmpy2.exp(1)

    elif seed_type == "misc":
        z = gmpy2.exp(gmpy2.sqrt(2))
    return(z)


#--- 5. Main

# set u=1, v=0 to show all k from k=n-start down to k=n-start-iters
# to show results only for k=v mod u, try u=60, v=25 
u = 1 ## 60 ## 60 # integer
v = 0 ## 25 # residue modulo u
# n = 3*7*11*13*u + v # choose n such that  n = u mod v 
n = 100000
truncate = True
start = 0
iters = 50000
iters = min(n-start, iters)
ctx = gmpy2.get_context()  
ctx.precision = 2*n 

seed_type = "primorial"

ncolors = 6  #  try number with many divisors: 12, 30, ...

colorTable = generate_contrasting_colors(ncolors)
z = initialize_reverse_seed(seed_type)
(arr_count1, arr_colors, xvalues) = iexp(n, start, iters, ncolors, z, u, v, truncate)


#--- 6. Create the main plot

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.facecolor'] = 'black'

plt.scatter(xvalues, arr_count1,s=0.2, c=arr_colors) 
#plt.plot(xvalues, arr_count1,linewidth=0.1, c='gold') 

plt.axhline(y=n/2,color='red',linestyle='--',linewidth=0.6,dashes=(5,10))   
plt.axhline(y=n/5, color='black', linestyle='--', linewidth = 0.6, dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))

for k in range(1,15):
    plt.axvline(x=k*n/(k+1),c='gray',linestyle='--',linewidth=0.6,dashes=(5, 10))

plt.ylim([-0.01 * n, 0.52*n]) 
plt.xlim([0.70 * n, 1.02*n]) 
plt.show()


#--- 7. Create the moving average plot

arr_avg = []
arr_xval = []
arr_count1 = np.array(arr_count1)
w = 60
for k in range(0, len(arr_count1)-w):
    #tmp = arr_count1[k:k+w]
    #print(k, len(tmp)) ##,arr_count1[k])
    y_avg = np.average(arr_count1[k:k+w])
    arr_avg.append(y_avg)
    arr_xval.append(n-start-k)

plt.scatter(arr_xval,arr_avg,s=0.02, c='gold')
plt.axhline(y=n/2,color='red',linestyle='--', linewidth=0.6,dashes=(5,10))   
plt.axhline(y=n/5, color='black', linestyle='--', linewidth = 0.6, dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))
plt.xlim([0.50*n, 1.01*n])
for k in range(1,15):
    plt.axvline(x=k*n/(k+1),c='gray',linestyle='--',linewidth=0.6,dashes=(5, 10))
plt.show() 

nv = len(arr_xval)
st = 0 ##int(4*n/5)
arr_delta =  np.array(arr_avg[0:nv-1])-np.array(arr_avg[1:nv])
plt.scatter(arr_xval[st+1:nv], arr_delta[st+0:nv-1], s=0.08, c=arr_colors[st+1:nv])
for k in range(1,15):
    plt.axvline(x=k*n/(k+1),c='gray',linestyle='--', linewidth=0.6,dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))
plt.axhline(y=0,color='red',linestyle='--', linewidth=0.6,dashes=(5,10))   
plt.show()


#--- 8. One more scatterplot

nv = len(arr_count1)
w = 1
arr_delta = np.array(arr_count1[0:nv-w])-np.array(arr_count1[w:nv])
plt.scatter(xvalues[w:nv], arr_delta, s=0.08, c=arr_colors[w:nv])
for k in range(1,15):
    plt.axvline(x=k*n/(k+1),c='gray',linestyle='--', linewidth=0.6,dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))
plt.axhline(y=0,color='red',linestyle='--', linewidth=0.6,dashes=(5,10))   
plt.show()
