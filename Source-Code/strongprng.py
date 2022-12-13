# By Vincent Granville, www.MLTechniques.com
 
import time
import random
import numpy as np

size =  400        # number of binary digits in each number
Niter =  5000      # number of quadratic irrationals
start = 0          # first value of (y, z) is (1, start)
yseed = 1          # y = yseed
offset =  100      #    skip first offset digits (all zeroes) of each number
PRNG = 'Quadratic' # options: 'Quadratic' or 'Mersenne' 
output = True      # True to print results (slow)

squareFreeList = {}
digits = {}
accepted = 0  # number of accepted seeds

for iter in range(start, Niter): 

    y = yseed # use variable  y
    z = iter
    c = (z - 1)**2 + 8*y 

    # represent c as a * b where a is square and b is square-free
    d = int(np.sqrt(c))
    a = 1 
    for h in range(2, d+1):
        if c % (h*h) == 0:  # c divisible by squared h
            a = h*h
    b = c // a   # integer division

    if b > 1 and b not in squareFreeList:
        q = (-(z - 1) + np.sqrt(c)) / 4  # number associated to seed (y, z); ~ y/(z-1)
        squareFreeList[b]=(y,z)          # accept the seed (y, z)
        accepted += 1

start = time.time()

for b in squareFreeList: 

    y = squareFreeList[b][0]
    z = squareFreeList[b][1]
         
    for k in range(size): 

        # trick to make computations faster
        y2 = y + y   
        y4 = y2 + y2
        z2 = z + z

        # actual computations
        if z < y2:
            y = y4 - z2
            z = z2 + 3
            digit = 1 
        else:
            y = y4
            z = z2 - 1
            digit = 0
        if k >= offset:
            digits[(b,k)] = digit

end = time.time()
print("Time elapsed:",end-start)

if output == True:
    OUT=open("strong4.txt","w")
    separator="\t"  # could be "\t" or "" or "," or " "
    if PRNG == 'Mersenne':
        random.seed(205) 
    for b in squareFreeList:
        OUT.write("["+str(b)+"]")
        for k in range(offset, size):
            key = (b, k)
            if PRNG == 'Quadratic':
                bit = digits[key] 
            elif PRNG == 'Mersenne':
                bit = int(2*random.random())
            OUT.write(separator+str(bit))
        OUT.write("\n")
    OUT.close()

print("Accepted seeds:",accepted," out of",Niter)
