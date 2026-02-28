import numpy as np
import gmpy2

def combinations_lexicographic_list(n, k):

    # Return a list with all k-combinations of range(n) in lexicographic order.
    # Each combination is a tuple of indices (0..n-1).

    if k < 0 or k > n:
        return []

    # initial combination: [0, 1, ..., k-1]
    c = list(range(k))
    cnt_00_min = 2*n
    cnt_00_max = -1
    cnt_01_min = 2*n
    cnt_01_max = -1
    cnt_10_min = 2*n
    cnt_10_max = -1
    cnt_11_min = 2*n
    cnt_11_max = -1
    flag = ''
    hash = {}
    print("\nFinding extremes, exhaustive search... \n")

    while True:

        str_0 = ''
        str_1 = ''
        for idx in range(n):
            if idx in c:
                str_0 += '1'
                str_1 += '0'
            else:
                str_0 += '0'
                str_1 += '1'

        x0 = two_third + gmpy2.mpz(str_0, 2)/2**n
        if x0 >= 1:
            x0 -= 1
        x0_bin = gmpy2.digits(x0, 2)[0]
        x0_bin = x0_bin[0:n]

        x1 = two_third + gmpy2.mpz(str_1, 2)/2**n
        if x1 >= 1:
            x1 -= 1
        x1_bin = gmpy2.digits(x1, 2)[0]
        x1_bin = x1_bin[0:n]

        cnt_00 = x0_bin.count('0')
        cnt_01 = x0_bin.count('1')
        cnt_10 = x1_bin.count('0')
        cnt_11 = x1_bin.count('1')
        
        combo = np.copy(c)

        if cnt_00 < cnt_00_min:
            cnt_00_min = cnt_00
            hash['min_00'] = (cnt_00, str_0, x0_bin, combo)
            flag += 'a'
        if cnt_00 > cnt_00_max:
            cnt_00_max = cnt_00
            hash['max_00'] = (cnt_00, str_0, x0_bin, combo)
            flag += 'b'
        if cnt_01 < cnt_01_min:
            cnt_01_min = cnt_01
            hash['min_01'] = (cnt_01, str_0, x0_bin, combo)
            flag += 'c'
        if cnt_01 > cnt_01_max:
            cnt_01_max = cnt_01
            hash['max_01'] = (cnt_01, str_0, x0_bin, combo)
            flag += 'd'

        if cnt_10 < cnt_10_min:
            cnt_10_min = cnt_10
            hash['min_10'] = (cnt_10, str_1, x1_bin, combo)
            flag += 'A'
        if cnt_10 > cnt_10_max:
            cnt_10_max = cnt_10
            hash['max_10'] = (cnt_10, str_1, x1_bin, combo)
            flag += 'B'
        if cnt_11 < cnt_11_min:
            cnt_11_min = cnt_11
            hash['min_11'] = (cnt_11, str_1, x1_bin, combo)
            flag += 'C'
        if cnt_11 > cnt_11_max:
            cnt_11_max = cnt_11
            hash['max_11'] = (cnt_11, str_1, x1_bin, combo)
            flag += 'D'

        if flag != '':
            print(c, str_0, x0_bin, cnt_00, cnt_01,
                        str_1, x1_bin, cnt_10, cnt_11, flag)
            flag = ''

        # find rightmost element that can be incremented
        i = k - 1
        while i >= 0 and c[i] == i + (n - k):
            i -= 1

        if i < 0:
            # all combinations generated
            break

        # increment this element
        c[i] += 1

        # reset the tail to the minimal increasing sequence
        for j in range(i + 1, k):
            c[j] = c[j - 1] + 1

    return(hash) 


#--- Main

n = 18 
k = 4
str = ''
for idx in range(n):
   if idx % 4 in (0,2):
       str += '1'
   else:
       str += '0'

ctx = gmpy2.get_context()  
ctx.precision = 2*n 
two_third = gmpy2.mpz(str, 2)/2**n
str2 = gmpy2.digits(two_third, 2)[0]
str2 = str2[0:n]
print("2/3, truncated string:",str2)

hash = combinations_lexicographic_list(n, k)
print("\nSummary\n")

for key in hash:

    value = hash[key]
    count = value[0]
    before = value[1]
    after = value[2]
    combo = value[3]
    combo = [f"{t}" for t in combo]
    combo = ' '.join(combo)
    rho = count/n
    print("%3d %3d | %s %s %s %3d %4.2f| %s" %(n, k, key, before, after, count, rho, combo))

