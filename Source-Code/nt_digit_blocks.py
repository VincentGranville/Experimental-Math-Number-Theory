import gmpy2 
import numpy as np

np.random.seed(454)

def compute_correct_digits(sum, pow2, n): 

    # to double-check our digits match those of theoretical limit
    bsum = bin(sum)[2:]
    bsum_u_exact = gmpy2.mpfr((pow2)/gmpy2.sqrt(pow2*pow2-4))
    bsum_main_exact = bsum_u_exact - bsum_u_exact/pow2
    bsum_constant_exact = pow2/gmpy2.mpfr(pow2*pow2-4)
    bsum_exact = bsum_main_exact + bsum_constant_exact
    bsum_exact = bin(gmpy2.mpz(pow2**(2*n) * bsum_exact))[2:]

    k = 0
    while k < len(bsum) and bsum[k] == bsum_exact[k]:
        k = k+1
    return(k) 


#--- 1. Main

ctx = gmpy2.get_context() 
n = 10000
nu = 1 
N = 100*n    # precision 
pow2 = 2**nu 
ctx.precision = N 
mode = 'SmallGrowth'  # options: 'Bernoulli', 'Fixed', 'SmallGrowth'

sum = 0
sum_y = 0
sum_u = 0
sum_v = 0
arr_q1 = []
arr_q2 = []
arr_q3 = []
arr_xval = []
arr_rlen = []
arr_delta2 = []
arr_delta3 = []

# sqrt2 is test number to compare digit distrib. with those of y and y'
sqrt2 = gmpy2.sqrt(2)  # this number used for comparison purpose
sqrt2  = gmpy2.mpz(2**(N) * sqrt2)  # turn sqrt2 into integer with N bits

for k in range(n):

    # B_max corresponds to M_{k, nu} in the paper

    if mode == 'Bernoulli':
        B_k = gmpy2.bincoef(2*k, k)
        B_max = 4**k  # known fact: B_k < B_max

    elif mode == 'Fixed':
        # Generate B_k with (nu + dx) bits, each bit is '1' with proba p
        # if dx > 0, there is carry over (more if dx large compared to nu)
        # dx can be < 0, but (nu + dx) most be >= 0
        # proportion of '1' in y is p only if dx=0 
        dx = 40   
        p = 0.25
        u = np.random.binomial(n=1, p=p, size=nu+dx)
        stri = [str(bit) for bit in u]
        stri = "".join(stri)
        B_k = int(stri, 2)
        B_max = 2**(nu + dx) - 1

    elif mode == 'SmallGrowth':
        cx = 3.5
        bx = 1
        B_k = int(bx * k**cx) 
        # try: B_k = 3**k >> k   # equal to int[(3/2)**k]
        B_max = 0
        if B_k > 0:
            bits = len(bin(B_k)) - 2
            B_max = 2**(bits+1)  # equal to 2^(int(log2 B_k))

    u_k = B_k
    v_k = B_max - B_k 
    s_k = pow2 * u_k + v_k 
    
    sum_u = pow2**2 * sum_u + u_k
    sum_v = pow2**2 * sum_v + v_k
    sum = pow2 * sum_u + sum_v     # for y'
    sum_y = pow2 * sum_y + B_k     # for y
    # equivalenly, sum = pow2**2 * sum + s_k 

    if k % 25 == 0:
        # computations needed to visualize intermediary results
        stri1 = bin(sum_y)[2:]  
        stri2 = bin(sqrt2)[2:]
        stri3 = bin(sum)[2:]
        rlen = -1
        correct_digits = -1
        if mode == 'Bernoulli' and nu > 1:
            # if correct_digits stalls, restart with larger N
            correct_digits = compute_correct_digits(sum, pow2, n)
            rlen = correct_digits/len(stri3)
            stri1 = stri1[0:correct_digits]
            stri2 = stri2[0:correct_digits]
            stri3 = stri3[0:correct_digits]
        else:
            stri2 = stri2[0:len(stri3)]
        q1 = stri1.count('1')/len(stri1)
        q2 = stri2.count('1')/len(stri2)
        q3 = stri3.count('1')/len(stri3)
        arr_q1.append(q1)
        arr_q2.append(q2)
        arr_q3.append(q3)
        arr_xval.append(k)
        arr_rlen.append(rlen)
        arr_delta2.append(abs(q2-0.5)*np.sqrt(k))
        arr_delta3.append(abs(q3-0.5)*np.sqrt(k))

        print("%5d %5d %6.4f %6.4f %6.4f %6.4f" % 
             (k, correct_digits, rlen, q1, q2, q3))  

print()
print(stri1[-60:])  # last 60 digits of y
print(stri3[-60:])  # last 60 digits of y'
print()

#-

# compute rescaled y and y'; factor is a power of nu
# works even when series diverges (Bernoulli with nu in {0, 1})
# if Bernoulli with nu = 0, then y' = 4/3 

y = 0
y_prime = 0
for k in range(60):
    y+= int(stri1[k])/2**k
    y_prime += int(stri3[k])/2**k
print("y  = %14.12f" % (y))
print("y' = %14.12f" % (y_prime))


#--- 2. Main plot: % of dgits equal to '1' as k increases

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['axes.facecolor'] = 'black'

plt.plot(arr_xval, arr_q1, linewidth = 0.4, c='lightblue')  # for y
plt.plot(arr_xval, arr_q2, linewidth = 0.8, c='green')  # for sqrt2
plt.plot(arr_xval, arr_q3, linewidth = 0.8, c='red')  # for y'
plt.axhline(y=0.50,color='gold',linestyle='--', linewidth=0.6,dashes=(5,10))   

legend = plt.legend(["y","sqrt2","y'"])
plt.setp(legend.get_texts(), color='white')
# plt.ylim(0.48, 0.52)
# plt.xlim(1000, n)
plt.show()


#--- 3. Other plots

if mode == 'Bernoulli':
    # show proportion of correct digits obtained after k iter
    plt.plot(arr_xval, arr_rlen, linewidth = 0.8)
    plt.show()

# show how far y' is to having 50% of `1' at iter k
plt.plot(arr_xval, arr_delta2, linewidth = 0.8, c='green')
plt.plot(arr_xval, arr_delta3, linewidth = 0.8, c='red')
plt.axhline(y=0.0,  color='gold',linestyle='--', linewidth=0.6,dashes=(5,10))   
legend = plt.legend(["sqrt2","y'"])
plt.setp(legend.get_texts(), color='white')
#plt.xlim(1000, n)
plt.show()
