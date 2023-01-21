import matplotlib.pyplot as plt
import mpmath
import numpy as np
from primePy import primes

# On WolframAlpha: DirichletL[4,2,s], s = sigma + it
#     gives the Dirichlet L-function  for the Dirichlet character  with modulus k and index j.
#
# References:
#     https://www.maths.nottingham.ac.uk/plp/pmzcw/download/fnt_chap4.pdf
#     https://mpmath.org/doc/current/functions/zeta.html
#     f(s) = dirichlet(s, [0, 1, 0, -1]) in MPmath
#     https://mathoverflow.net/questions/439018/zeros-of-dirichlet-function-ls-chi-4/439019#439019

m = 150000 # max number of primes to include in the Euler prduct
p1 = []
p3 = []
p  = []
cnt1 = 0
cnt3 = 0
cnt  = 0
for k in range(m):
    if primes.check(k) and k>1:
        if k % 4 == 1:
            p1.append(k)
            p.append(k)
            cnt1 += 1
            cnt += 1
        elif k % 4 ==3:
            p3.append(k)
            p.append(k)
            cnt3 += 1
            cnt += 1

cnt1 = len(p1)
cnt3 = len(p3)
n = min(cnt1, cnt3)
max = min(p1[n-1],p3[n-1])

print(n,p1[n-1],p3[n-1])
print()

sigma = 0.95
t_0 = 6.0209489046975965 # 0.5 + t_0*i is a root of DL4

DL4 = []
imag = []
print("------ MPmath library")
for t in np.arange(0,1,0.25):
    f = mpmath.dirichlet(complex(sigma,t), [0, 1, 0, -1]) 
    DL4.append(f)
    imag.append(t)
    r = np.sqrt(f.real**2 + f.imag**2)
    print("%8.5f %8.5f %8.5f" % (t,f.real,f.imag))

print("------ scrambled product")  
for t in np.arange(0,1,0.25):
    prod = 1.0
    for k in range(n):
        num1 = 1 - mpmath.power(1/p1[k],complex(sigma,t))
        num3 = 1 + mpmath.power(1/p3[k],complex(sigma,t))
        prod *= (num1 * num3)
    prod = 1/prod
    print("%8.5f %8.5f %8.5f" % (t,prod.real,prod.imag))

DL4_bis = []
print("------ scrambled swapped") 
for t in np.arange(0,1,0.25):
    prod = 1.0
    for k in range(n):
        num1 = 1 + mpmath.power(1/p1[k],complex(sigma,t))
        num3 = 1 - mpmath.power(1/p3[k],complex(sigma,t))
        prod *= (num1 * num3)
    prod = 1/prod
    DL4_bis.append(prod)
    print("%8.5f %8.5f %8.5f" % (t,prod.real,prod.imag))

#plt.plot(x,y,c='blue')
#plt.plot(uu,vv,c='blue')
#plt.show()

print("------ compare zeta with DL4 * DL4_bis")
for i in range(len(DL4)):
    t = imag[i]
    if t == 0 and sigma == 0.5:
        print("%8.5f" % (t))
    else:
        zeta = mpmath.zeta(complex(2*sigma,2*t))
        prod = DL4[i] * DL4_bis[i] / (1 - 2**(-complex(2*sigma,2*t))) 
        print("%8.5f %8.5f %8.5f %8.5f %8.5f" % (t,zeta.real,zeta.imag,prod.real,prod.imag))

print("------ correct product")
for t in np.arange(0,1,0.25):
    prod = 1.0
    chi = 0
    k = 0
    while p[k] <= max:
        pp = p[k]
        if pp % 4 == 1:
            chi = 1
        elif pp % 4 == 3:
            chi = -1
        num = 1 - chi * mpmath.power(1/pp,complex(sigma,t))
        prod *= num
        k = k+1
    prod = 1/prod
    print("%8.5f %8.5f %8.5f" % (t,prod.real,prod.imag))

print("------ series")
for t in np.arange(0,1,0.25):
    sum = 0.0
    flag = 1
    k = 0
    while 2*k + 1 <= 10000:
        num = flag * mpmath.power(1/(2*k+1),complex(sigma,t))
        sum = sum + num
        flag = -flag
        k = k + 1
    print("%8.5f %8.5f %8.5f" % (t,sum.real,sum.imag))
