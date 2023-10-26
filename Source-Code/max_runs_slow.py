import math
import gmpy2

# requirement: 0 < p < q
p = 1 
q = 2 
x0 = math.sqrt(p/q)
N = 1000 # precision, in number of binary digits for x0 

# compute and store in bsqrt (a string) the N first binary digits of x0 = sqrt(p/q)
base = 2
bsqrt = gmpy2.isqrt( (2**(2*N) * p) // q ).digits(base) 

for n in range(1, N):

    if n == 1:
        u = p * 4**n  
        v = int(x0 * 4**n) 
        if v % 2 == 0:
            v = v - 1
    else: 
        u = 4*u
        v = 2*v + 1
    steps = 0
    while q*v*v < u:
        v = v + 2
        steps += 1   # steps is always 0, 1, or 2
    v = v - 2
    delta = u - q*v*v 
    d = bsqrt[n-1]    # binary digit of x0 = sqrt(p/q), in position n 
    
    ## delta2 = delta >> (n - 1)
    ## res = 5/2 + n - math.log(delta,2) - math.log(n, 2)

    run = int(n + 1 + math.log(p*q, 2)/2 - math.log(delta, 2) ) 
    if d == "0" or run == 0:
        run = "" 

    print("%6d %1s %2s %1d" % (n, d, str(run), steps))
