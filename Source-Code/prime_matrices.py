import numpy as np

#--- 1. For each determinant det, create list of matrices and count them

hash_det = {}
hash_det_list = {}
max_det = 1000
upper_L = 2+int(max_det/3)

for c in range(0, upper_L):           
    if c % 50 == 0:
        print("c = %4d / %4d" % (c, upper_L))
    for b in range(c, min(upper_L,max_det-c)):
        sum = b**3 + c**3 
        prod = -3*c*b
        for a in range(b, min(upper_L, max_det-b-c)): 
            det = sum + a*(a*a + prod)  
            if 0 < det < max_det:
                if det in hash_det:
                    hash_det[det] += 1
                    dlist = hash_det_list[det] 
                    dlist = (*dlist, [a, b, c])
                    hash_det_list[det] = dlist
                else:
                    hash_det[det] = 1
                    hash_det_list[det] = ([a, b, c],)

mode = 'show_all' 
# options: 'show_nontrivial_matrices', 'show_all'

for det in range(1, max_det):
    if det not in hash_det:
        hash_det[det] = 0
        hash_det_list[det] = ()

def frequency_cumul(det_count):

    # For each Delta with |S_Delta| == det_count, show all matrices with determinant Delta 
    # |S_Delta| is the number of elements in S_Delta; Delta is denoted as det in the code
    # Also gather cumulative stats to further show plots 

    arr_det = []
    arr_count = []
    sum = 0

    for det in range(1, max_det):
        if hash_det[det] == det_count:
            arr_det.append(det)
            sum += 1 
            arr_count.append(sum)
            alist = hash_det_list[det]
            if mode == 'show_all':
                S_Delta = alist
            else:
                S_Delta = ()
                for Matrix in alist:
                    llambda = Matrix[0] + Matrix[1] + Matrix[2]
                    if llambda != det: 
                        S_Delta = (*S_Delta, Matrix)
            if len(S_Delta) > 0:
                print(det, hash_det[det], S_Delta)
        else:
            arr_det.append(det)
            arr_count.append(sum)

    return(arr_det, arr_count)

#--- 2. plot results, print matrix list and count for each determinant

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 'x-small' 

print()
for det_count in range(0,7):
    arr_det, arr_count = frequency_cumul(det_count)
    plt.plot(arr_det, arr_count, linewidth = 0.5, label = r'k =%d' %det_count)
plt.legend(loc="upper left")
plt.show()


#--- 3. Compute products

def amatrix(a, b, c):
    M = [[a, b, c],
         [c, a, b],
         [b, c, a]]
    M = np.array(M)
    return(M)

# lists of matrices with given determinant
A_det = 108
B_det = 189
Alist = hash_det_list[A_det] 
Blist = hash_det_list[B_det] 
print()

# compute all cross-products from Alist with Blist
canonical_form = False
for A in Alist:
    for B in Blist:
        if canonical_form:
            A = sorted(A, reverse=True)
            B = sorted(B, reverse=True)
        M = amatrix(A[0], A[1], A[2]) @ amatrix(B[0], B[1], B[2])
        M = M[0] 
        if canonical_form:
            M = sorted(M, reverse=True)
        M = np.array(M)
        print(A,"*", B, "=", M)

