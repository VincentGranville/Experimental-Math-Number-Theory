# Faster version than number_theory_fast_v2.py
#   - at iteration k, keep only 2n-k digits in S(n, k, x) instead of 2n
#   - also remove the trailing 0 on the right, in S(n, k, x)
#   - drawback: I get 19985 correct digits instead of 19998 if n = 20000

n = 20000 
H = int(1.15*n)  

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

# precision set to L bits to keep at least about n correct bits till k=kmax
ctx = gmpy2.get_context()  
L = n + int(kmax+1) 
ctx.precision = L     

# p = +1: leading to e (Euler's number), about n correct bits after n iter
#         seed with n+1 bits, all 0 except rightmost and leftmost
# p = +3: leading to cube of e, about n correct bits after n iter
#         seed with n+1 bits, all 0 except 2 rightmost and leftmost
# p = -1: leading to inverse of e, about n correct bits after n iter
#         seed with n bits, all 1
# p must be an integer != 0; use p=0 for random seed

p = 1   # try -1, 0 (random seed), 1, 3 

def create_random_seed(n, cnt1):

    # create random seed of length n-1 with cnt1 '1' at random locations
    # and add a 1' at both ends; cnt1 must be <= n-1
 
    cnt1 = min(cnt1, n-1)
    numpy_seed = 453    # numpy seed to initiate numpy PRNG, not the model seed
    np.random.seed(numpy_seed)
    random_locations = np.random.choice(np.arange(1, n),size=cnt1,replace=False)
    prod = 2**n + 1     # seed with n+1 '0' except a '1' at both ends
    for position in random_locations:
        # add the cnt1 random '1's between both ends
        prod += 2**int(position)
    return(prod)


# create seed with n+1 bits if p>=0, or n bits if p<0
if p != 0:
    prod = gmpy2.mpz(2**n + p)           
else:
    # random seed with number of '1' in seed to cnt1+2
    # test: set n = 30000; cnt1 = 0, 3, 4, 5, 100  and see what happens!
    cnt1 = 3  
    prod = create_random_seed(n, cnt1)  

# local variables
arr_count1 = []
arr_colors = []
xvalues = []
ecnt1 = -1

OUT = open("digit_sum.txt", "w")

for k in range(1, H+1): 

    prod = prod*prod
    pstri = bin(int(prod))
    stri = pstri[0:2*n-k]      ### faster than pstri[0: L+2] in older version
    stri = rstrip_zeros(stri)  ### new to this version (faster)
    prod = int(stri, 2)
    prod = gmpy2.mpz(prod) 

    if k > kmin and k < kmax:   
        stri = stri[2:]
        if k == n:
            e_approx = stri
        estri = stri[0:n]   # leftmost n digits
        ecnt1 = estri.count('1')    ### old version
        ### ecnt1 = stri.count('1') ### new version
        arr_count1.append(ecnt1)
        color = assign_color(k)
        arr_colors.append(color)
        xvalues.append(k)
        OUT.write(str(k) + "\t" + str(ecnt1) + "\n")
 
    if k%1000 == 0:
        print("%3d %3d" %(k, ecnt1))
    if stri[-1] == '0':
        print(k, stri[-10:])

OUT.close()


#--- 2. Compute bits of e and count correct bits in my computation

from mpmath import mp

# Set precision to L binary digits
mp.dps = int(L*np.log2(10))
e_value = (mp.e)**abs(p)  # Get e^|p| in decimal

if p > 0:
    # Convert e_value to binary string
    e_binary = bin(int(e_value * (2 ** n)))[2:] 

elif p < 0: 
    e_iapprox = int(e_approx, 2)  # convert string e_approx to integer
    e_ivalue = int(2**(2*n) * e_value)
    one = e_iapprox * e_ivalue
    e_approx = bin(one)[2:]  
    e_binary = "1" * (2*n)  # string of 2n bits, all '1'

if p != 0:
    k = 0
    while e_approx[k] == e_binary[k]:
        k += 1
    # e_binary should be equal to e_approx up to about n bits  
    print("\n%d correct digits (n = %d)" %(k, n))


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

if p > 0:
    # we start with about 0% of 1 going up to about 50%
    ymax = 0.52 * n
    plt.ylim([-0.01 * n, ymax]) 
elif p < 0:
    # we start with 100% of 1 going down to about 50%
    ymax = 1.01 * n
    plt.ylim([0.40 * n, ymax])
elif p == 0: 
    ymax = 1.00 * n
    plt.ylim([-0.01 * n, ymax])
##plt.xlim([kmin, kmax])
plt.xlim([0.7*n, kmax])


plt.show()

#--- 4. Create the average plot

arr_avg = []
arr_xval = []
arr_count1 = np.array(arr_count1)

for k in range(0, int(kmax-12)):
    y_avg = np.average(arr_count1[k:k+12])
    arr_avg.append(y_avg)
    arr_xval.append(k)

plt.scatter(arr_xval,arr_avg,s=0.0002, c='gold')
plt.axhline(y=n/2,color='red',linestyle='--',linewidth=0.6,dashes=(5,10))   
plt.axhline(y=n/5, color='black', linestyle='--', linewidth = 0.6, dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.6, dashes=(5, 10))
plt.xlim(0.00*n, kmax-12)
plt.ylim(-0.01*n, ymax)

for k in range(1,15):
    plt.axvline(x=k*n/(k+1),c='gray',linestyle='--',linewidth=0.6,dashes=(5, 10))

plt.show() 
