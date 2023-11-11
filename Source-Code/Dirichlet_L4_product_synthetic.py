import numpy as np

N = 10000000
p = 1.00 
q = 0.90  
stdev = 0.50
seed = 564
np.random.seed(seed)
start = 20

u = 0
v = 0
a = np.zeros(N)
b = np.zeros(N)

for n in range(2, N):

    u += -0.5 + np.random.randint(0, 2) 
    v += np.random.normal(0, stdev)/n**q  
    a[n] = u / n**p 
    b[n] = v

    if n % 50000 == 0:
        sa = np.std(a[start:n])
        sb = np.std(b[start:n])
        r = sa / sb
        c = r * b[n] - a[n]
        print("n = %7d r =%8.5f an =%8.5f bn =%8.5f c =%8.5f sa =%8.5f sb=%8.5f" 
                  %(n, r, a[n], b[n], c, sa, sb)) 
    
