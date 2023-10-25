import gmpy2

p = 1 
q = 2 
N = 1000000000 # precision, in number of binary digits 

# compute and store in bsqrt (a string) the N first binary digits of sqrt(p/q)
base = 2
bsqrt = gmpy2.isqrt( (2**(2*N) * p) // q ).digits(base) 

last_digit = -1
L = 0
max_run = 0

for n in range(0, N):
    d = int(bsqrt[n])  # binary digit
    if d == 0:
        L += 1
    if d == 1 and last_digit == 0:
            run = L 
            if run > max_run:
                max_run = run
                print(n-L, run, max_run)
            L = 0
    last_digit = d
