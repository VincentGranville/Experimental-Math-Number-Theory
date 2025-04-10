# Faster version than number_theory_fast_v2.py
#   - at iteration k, keep only 2n-k digits in S(n, k, x) instead of 2n
#   - also remove the trainling 0 on the right, in S(n, k, x)
#   - drawback: I get 19985 correct digits instead of 19998 if n = 20000

from primePy import primes
import gmpy2
import numpy as np

n = 100000 
H = int(1.1*n)  

import colorsys 

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

ncolors = 6  #  try number with many divisors: 12, 30, ...
colorTable = generate_contrasting_colors(ncolors)


def rstrip_zeros(string):

    # remove '0' on the right after last '1'

    newstring = string
    if string[-1] == '0':
        k = -1
        while string[k] == '0':
            k -= 1
        newstring = string[:k+1]
    return(newstring)


#--- 1. Main

import gmpy2
import numpy as np

kmin = 0.00 * n  # do not compute digit count if k <= kmin
kmax = 1.15 * n  # do not compute digit count if k >= kmax
kmax = min(H, kmax)

# precision set to L bits to keep at least about n correct bits till k=n
ctx = gmpy2.get_context()  
ctx.precision = 2*n 

# p = 2*3*5*7*11*13*17*19*23*29
p = 1   # denoted as x in the paper

# first n binary digits at iteration k=n are those of exp(p) 
# if p irrational, seed = 2^(2n) + int(2^n * p)
# if p integer, seed = 2^n + p

iplog = 0

if p != int(p):
    # for p irrational, like p = sqrt(2)
    p = gmpy2.floor((2**n)*gmpy2.mpfr(p)) 
    prod = gmpy2.floor(2**(2*n) + p)
else:
    # for integer, small or large
    iplog = gmpy2.floor(gmpy2.log2(abs(p)))
    prod = gmpy2.floor(2**n + p)  

# local variables
arr_count1 = []
arr_colors = []
xvalues = []
ecnt1 = -1
e_approx = "N/A"

OUT = open("digit_sum.txt", "w")

for k in range(1, H+1): 

    prod = prod*prod
    pstri = bin(gmpy2.mpz(prod))  # mpz is round to integer, not floor
    stri = pstri[0:2*n-k]         # faster than pstri[0: L+2] in older version
    stri = rstrip_zeros(stri)     # new to this version (faster)
    prod = int(stri, 2)
    prod = gmpy2.floor(prod) 

    if k > kmin and k < kmax:   
        stri = stri[2:]
        lstri = len(stri)
        if k == n-iplog:
            e_approx = stri
        estri = stri[0:n]   # leftmost n digits
        ecnt1 = estri.count('1')   
        ecnt1f = stri.count('1') 
        arr_count1.append(ecnt1)
        color = colorTable[k % ncolors] 

        arr_colors.append(color)
        xvalues.append(k)
        OUT.write(str(k)+"\t"+str(ecnt1)+"\t"+str(lstri)+"\t"+str(ecnt1f)+"\n")
 
        if k%1000 == 0:
            print("%6d %6d %6d %6d" %(k, ecnt1, lstri,ecnt1f))
    if stri[-1] == '0':
        print(k, stri[-10:])

OUT.close()


#--- 2. Compute bits of e and count correct bits in my computation

# Set precision to L binary digits
gmpy2.get_context().precision = 4*n
if p == int(p):  
    e_value = gmpy2.exp(p/2**iplog) 
else:
    e_value = gmpy2.exp(p)

# Convert e_value to binary string
e_binary = gmpy2.digits(e_value, 2)[0]

k = 0
while e_approx[k] == e_binary[k]:
        k += 1
# e_binary should be equal to e_approx up to about n bits 
if p == int(p): 
    print("\n%d correct digits (n = %d, iplog = %d)" %(k, n, iplog))

e_approx_decimal = 0
for k in range(80): 
    e_approx_decimal += int(e_approx[k])/(2**k) 
print("e_approx, up to power of 2:", e_approx_decimal)


#--- 3. Create the main plot

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.facecolor'] = 'black'

plt.scatter(xvalues, arr_count1,s=0.01, c=arr_colors) 
#plt.plot(xvalues, arr_count1,linewidth=0.1, c='gold') 

plt.axhline(y=n/2,color='red',linestyle='--',linewidth=0.6,dashes=(5,10))   
plt.axhline(y=n/5, color='black', linestyle='--', linewidth = 0.6, dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))

for k in range(1,15):
    plt.axvline(x=k*n/(k+1),c='gray',linestyle='--',linewidth=0.6,dashes=(5, 10))

# we start with about 0% of 1 going up to about 50%
ymax = 0.52 * n
plt.ylim([-0.01 * n, ymax]) 
#plt.ylim([-0.01 * n, 1.01*n]) 

plt.xlim([kmin, kmax]) 
# plt.xlim([0.0*n, kmax])
plt.show()


#--- 4. Create the moving average plot

arr_avg = []
arr_xval = []
arr_count1 = np.array(arr_count1)
w = 6 

for k in range(0, len(arr_count1)-w):
    y_avg = np.average(arr_count1[k:k+w])
    arr_avg.append(y_avg)
    arr_xval.append(k)

plt.scatter(arr_xval,arr_avg,s=0.0002, c='gold')
plt.axhline(y=n/2,color='red',linestyle='--',linewidth=0.6,dashes=(5,10))   
plt.axhline(y=n/5, color='black', linestyle='--', linewidth = 0.6, dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))
plt.xlim(0.00*n, kmax-w)
plt.ylim(-0.01*n, ymax)

for k in range(1,15):
    plt.axvline(x=k*n/(k+1),c='gray',linestyle='--',linewidth=0.6,dashes=(5, 10))

plt.show() 

nv = len(arr_xval)
st = int(4*n/5)
arr_delta = np.array(arr_avg[1:nv]) - np.array(arr_avg[0:nv-1])
plt.scatter(arr_xval[st+1:nv], arr_delta[st+0:nv-1], s=0.08, c=arr_colors[st+1:nv])
for k in range(1,15):
    plt.axvline(x=k*n/(k+1),c='gray',linestyle='--',linewidth=0.6,dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))
plt.axhline(y=0,color='red',linestyle='--',linewidth=0.6,dashes=(5,10))   
plt.xlim(st, nv)
plt.show()


#--- 5. Create AR scatterplot

nv = n ###  len(arr_xval)
plt.scatter(arr_count1[nv-2000-1:nv-1],arr_count1[nv-2000:nv],s=0.04, c=arr_colors[nv-2000-1:nv-1])
plt.show()

