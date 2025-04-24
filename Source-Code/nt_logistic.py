import gmpy2
import numpy as np

n = 99999 # choose n divisible by 3 if x = 9 and ncolors = 6
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

ncolors = 6 # 4 colors for hybrid case, 6 for quantic
colorTable = generate_contrasting_colors(ncolors)


#--- 1. Main

import gmpy2
import numpy as np

kmin = 0.00 * n  # do not compute digit count if k <= kmin
kmax = 1.15 * n  # do not compute digit count if k >= kmax
kmax = min(H, kmax)

# precision set to L bits to keep at least about n correct bits till k=n
ctx = gmpy2.get_context()  
ctx.precision = 2*n 

# x = gmpy2.mpz(2*3*5*7*11*13*19*23*29) 
x = gmpy2.mpz(1) 

prod = gmpy2.mpfr(x/2**(2*n))

# local variables
arr_count1 = []
arr_colors = []
xvalues = []
ecnt1 = -1
e_approx = "N/A"

OUT = open("digit_sum.txt", "w")

for k in range(1, H+1): 

    prod = 4*prod*(1 - prod)
    pow = 2*n  # 2n --> 0.1111.. | 4n --> 1.0000.. 
    pstri = bin(gmpy2.mpz(2**(pow) * prod)) 
    stri = pstri[0:2*n]

    if k > kmin and k < kmax:   
        stri = stri[2:]
        lstri = len(stri)
        if k == n:
            e_approx = stri
        # estri = stri[0:n]             # for digit sum 
        estri = stri[max(0,n-k):2*n-k]  # for adjusted digit sum
        ecnt1 = estri.count('1') * n / (1+len(estri))   
        arr_count1.append(ecnt1)
        color = colorTable[k % ncolors] 

        arr_colors.append(color)
        xvalues.append(k)
        OUT.write(str(k)+"\t"+str(ecnt1)+"\t"+str(lstri)+"\n")
 
        if k%1000 == 0:
            print("%6d %6d %6d" %(k, ecnt1, lstri))

OUT.close()


#--- 2. Compute bits of sin^2(sqrt(x_) and count correct bits in my computation

# Set precision to L binary digits
gmpy2.get_context().precision = 4*n
e_value = gmpy2.sin(gmpy2.sqrt(x)) 
e_value = e_value * e_value

# Convert e_value to binary string
e_binary = gmpy2.digits(e_value, 2)[0]

k = 0
while k < len(e_approx) and e_approx[k] == e_binary[k]:
        k += 1
# e_binary should be equal to e_approx up to about n bits 
print("\n%d correct digits (n = %d)" %(k, n))

e_approx_decimal = 0
for k in range(0,80): 
    e_approx_decimal += int(e_approx[k])/(2**(k+1)) 
print("e_exact : %16.14f" % (e_value))
print("e_approx: %16.14f" % (e_approx_decimal))
print("Up to factor 2 at integer power of 2.")

#--- 3. Create the main plot

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.facecolor'] = 'black'

plt.scatter(xvalues, arr_count1,s=0.005, c=arr_colors) 
plt.axhline(y=n/2,color='red',linestyle='--', linewidth=0.6,dashes=(5,10))   
plt.axhline(y=n/5, color='black', linestyle='--', linewidth = 0.6, dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))

for k in range(1,15):
    plt.axvline(x=k*n/(k+1),c='gray',linestyle='--', linewidth=0.6,dashes=(5, 10))

# we start with about 0% of 1 going up to about 50%
plt.ylim([0.44*n, 1.01*n]) 
plt.xlim([0.45*n, 1.02*n])
plt.show()


#--- 4. Create AR scatterplot

nv = n
lag = 12 
tail = 2000
plt.scatter(arr_count1[n-tail-lag:n-lag], arr_count1[n-tail:n], s=0.4, c=arr_colors[n-tail-lag:n-lag])
# plt.plot(arr_count1[n-tail-lag:n-lag], arr_count1[n-tail:n], linewidth=0.6) ##, c=arr_colors[n-tail-lag:n-lag])
plt.axhline(y=n/2,color='red',linestyle='--',linewidth=0.6,dashes=(5,10))   
plt.axvline(x=n/2, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))

plt.show()

