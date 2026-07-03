# import sys
# sys.dont_write_bytecode = True  

import importlib
import gilbreath_lib as gil
importlib.reload(gil)

import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 'x-small'


#--- 1. Showing how central function gil.check works

nprimes = 30 
from primePy import primes
prime_list = primes.first(nprimes)
arr_a = prime_list[0:nprimes]
# arr_a = [2, 3, 5, 11, 13, 15, 31, 33, 35, 39, 51 ]
print("--- Part 1: gil.check\n")

flag, first_amax_idx, hash_stats = gil.check(arr_a, show_triangle = 'Full', mode = 'Full', delta = 'standard') #######################

right_diagonal = hash_stats['right_diagonal']
left_diagonal = hash_stats['left_diagonal']
failing_level = hash_stats['failing_level']
print("Right diagonal:", right_diagonal)
print("Left diagonal:", left_diagonal)
print("Failing_level:", failing_level)
print()


#--- 2. Detect failed sequences by simulation [Poisson increments]

parameters = { 
               'prime': [1.0, 1.3, 1.50, 5],
               'power': [0.9, 1.1, 1.25, 5, 1.50],
             }

model = 'prime'

n = 60  # number of values to add after arr_init
nsamples = 400 
llambda = 1.5 
rejection_sampling = True   # set to True to guarantee sequences are in the corridor
arr_init = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]  
rng = np.random.default_rng(seed=42) 

offset = len(arr_init)

from primePy import primes
prime_list = primes.first(n + offset)
p_n = prime_list[-1]
success = 0
in_corridor = 0
sum = 0
hash_fail = {}
print("--- Part 2: Simulations\n")

for sample in range(nsamples):

    arr_a = gil.generate_random_sequence(n, arr_init, llambda, parameters, model, rng)
    sum += arr_a[-1]
    flag, first_amax_idx, hash_stats = gil.check(arr_a, show_triangle = 'None')
    failing_level = hash_stats['failing_level']
    first_diff = hash_stats['first_diff']  
    left_diaginal = hash_stats['left_diagonal']
    valid, bounded, smooth, fvals, sequence = gil.admissibility(first_diff, parameters, model)  # fvals  empty if in corridor

    if flag == 'Fail' and (bounded or smooth):

        snapshot = arr_a[offset-1 : failing_level]
        left_diagonal = hash_stats['left_diagonal']
        sigma = 1 + left_diagonal[-1]
        failing_diff = tuple(first_diff[0:failing_level])
        gamma = np.max(failing_diff)
        hash_fail[failing_diff] = (failing_level, sigma, gamma, failing_diff[-2], failing_diff[-1])
        print("Failed:", sample, failing_level, sigma, bounded, smooth, arr_a[-1], 
                         p_n, failing_level, first_diff[0:failing_level]) 

    if  bounded and smooth:
        in_corridor += 1
    if flag == 'Success':
        success += 1
    if sample % 500 == 0:
        print("Sample", sample, success-1)

avg_last_value = sum/nsamples

for first_diff in hash_fail:

    value = gil.make_readable(hash_fail[first_diff])
    first_diff = gil.make_readable(first_diff)
    first_diff_seq = first_diff[offset-6:len(first_diff)]
    banned_found, banned_pattern = gil.banned_found_in_seq(first_diff_seq)
    banned_pattern = gil.make_readable(banned_pattern)

failed_deduped = len(hash_fail)
print("Simulations:", success, "/", nsamples, avg_last_value, p_n, in_corridor, failed_deduped)
print()

#--- 3. Canonical form of sequence 

n_levels = 400 
fail_levels = { 120:3, 281:3 }  # err1.png 
fail_levels = { 120:3, 143:3 }  # err2.png
fail_levels = { 120:3, 144:3 }  # err3.png

N = 1000 
primes_ = 'real'  # options: 'real' or 'fake'

if primes_ == 'real':
    prime_list = primes.first(N)
else:
    arr_init = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]  
    llambda = 2.5 
    prime_list = gil.generate_random_sequence(N, arr_init, llambda, parameters, model, rng)

s1 = np.copy(prime_list)
valid = False
print("--- Part 3: Canonical form\n")

s2_failed, s2_failed_gaps, s1_bottom, s1_max_gap, s1_proba = gil.canonical(s1, n_levels, fail_level = fail_levels, valid = valid)
s2_succes, s2_succes_gaps, s1_bottom, s1_max_gap, s1_proba = gil.canonical(s1, n_levels, valid = valid)

plt.gca().set_facecolor('black')
arr_k = np.arange(2, len(s1))
plt.plot(arr_k, np.log(s2_failed[2:]), c='red', lw=0.9)
plt.plot(arr_k, np.log(s2_succes[2:]), c='green', lw=0.6) 
for k in fail_levels:
    plt.axvline(x=k, color='gray', linestyle='--', linewidth=0.4)
plt.axvline(x=n_levels, color='gray', linestyle='--', linewidth=0.4)
plt.show()

print("Summary:", N, n_levels, np.max(s2_succes_gaps), np.max(s1_bottom), np.max(s2_succes), np.max(s1), s1_max_gap)
print("Probas:", s1_proba[0:5])

 
#--- 4. Find bounds for next prime, to keep the sequence successful moving forward

sequence_type = 'prime_gaps' # options: 'primes', 'log_gaps', 'power_gaps'
prime_type = 'real'  # 'real', 'simulated'
nprimes = 5000 
start = 2 # must be 2 or higher

test_mode = False  # very slow if True, will test conjecture
if test_mode:
    conjecture = 'satisfied' # satisfied until proven wrong
else:
    conjecture = 'untested'

if prime_type == 'real':
    from primePy import primes
    prime_list = primes.first(nprimes)
else:
    arr_init = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] 
    llambda = 2.5 # 5
    prime_list = gil.generate_random_sequence(nprimes, arr_init, llambda, parameters, model, rng)
    flag, acheck, hash_stats = gil.check(prime_list, show_triangle = 'None', mode = 'Fast')
    print("check status of generated sequence:", flag)

arr_k = []
arr_min = []
arr_max = []
arr_a = []
arr_old_a = []
resets_list = ()
a = prime_list[start-1]
print("\n--- Part 4: find q_plus, q_minus\n")

OUT = open("gil_output.txt", "wt") 
print("idx", "reset", "rvbool", "k", "a", "a-old_a", "max_a-a-k", "end", "counts", 
      "dmax", "dmax2", "dmax2_idx", "extract", file = OUT, sep="\t")

# mode = 'Full' to not stop as soon as success is detected and get full diagonal
flag, acheck, hash_stats = gil.check(prime_list[0:start], mode = 'Full')
right_diagonal = hash_stats['right_diagonal'] 
right_diagonal = gil.make_readable(right_diagonal) 
hash_output = gil.collect_statistics(right_diagonal)

start_time = time.perf_counter()

for k in range(start, nprimes): 

    # given previous primes, find possible range for next one to guarantee success
    # augmenting sequence with previous prime as starting point, seems to always succeeds 

    old_a = a  
    old_hash_output = hash_output
    
    max_a, trials1 = gil.max_gap_allowed(right_diagonal, old_a) 
    min_a, trials2 = gil.min_gap_allowed(right_diagonal, old_a)

    if sequence_type == 'prime_gaps': 
            a = prime_list[k]
    else:
        if sequence_type == 'log_gaps':
            # mimic prime growth 
            gap = int(np.log(max_a)) # max_a can be as large as 2a
        elif sequence_type == 'power_gaps':
            gap = int(max_a**0.35) # exponent must be < 1
        gap = 2 + (gap - gap % 2)
        a += gap

    arr_k.append(k)
    arr_min.append(min_a)
    arr_max.append(max_a)
    arr_a.append(a)
    arr_old_a.append(old_a)
        
    if test_mode:
        # check if all 'a' are in [old_a, max_a] succeed, and max_a + 2 fails
        satisfied, flag = gil.test_conjecture(old_a, max_a, right_diagonal)
        if not satisfied:
            conjecture = 'not satisfied'

    old_idx = old_hash_output['before_02_index']
    old_nu_2 = old_hash_output['nu_2']
    success, right_diagonal = gil.augment(right_diagonal, a)
    right_diagonal = gil.make_readable(right_diagonal) 
    hash_output = gil.collect_statistics(right_diagonal)

    reset     = hash_output['reset'] 
    counts    = hash_output['counts'] 
    dmax2     = hash_output['max_value'] 
    dmax2_idx = hash_output['max_index'] 
    end       = hash_output['before_02_value'] 
    idx       = hash_output['before_02_index'] 
    dmax      = hash_output['max_tail_value'] 
    # ri        = hash_output['reversal_idx'] 
    # rval      = hash_output['reversal_value']
    rvbool    = hash_output['reversal']

    if reset:
        resets_list = (*resets_list, k) 
    if k % 1000 == 0:
        print("progress:", k, "/", nprimes)

    extract = right_diagonal[1:20] 
    print(old_idx, idx, old_nu_2, reset, rvbool, k, a, a-old_a, max_a-a-k, end, counts,
           dmax, dmax2, dmax2_idx, extract, file = OUT, sep="\t")
    
end_time = time.perf_counter()
elapsed_ms = end_time - start_time
print(f"Elapsed time: {elapsed_ms:.2f} ms")
OUT.close()

gil. arr_plots(arr_k, arr_min, arr_max, arr_a) # , resets_list) 

# future project:  analyse p_n - 2 p_{n-1} + p_{n-2}


#--- 5. Produce trimmed sequences with Sieve of Eratosthenes

all_primes = False
max_period = 20000
if all_primes:
    from primePy import primes
    # use the first 100k primes
    moduli = primes.first(100000)
else:
    moduli = (2, 3, 5, 7, 11, 13, 17, 19)

period = 1
for p in moduli:
    period *= p - 1
print("period = ", period)
arr_a = [2, 3]
start = 4 
end =  min(period + len(moduli), max_period)

count = 0
k = start

while count < min(end, max_period): 
   keep = True
   index = 0
   while keep and index < len(moduli):
       mod = moduli[index]
       if k % mod == 0:
           keep = False
       index += 1
   if keep or k in moduli: 
       arr_a.append(k)
       count += 1
   k += 1
   if k % 1000000 == 0:
       print("/...", k, count)

print("\n--- Part 5: Sieve of Eratosthenes\n")
show_triangle = 'Compact'
arr_a_truncated = arr_a[0:max_period]
flag, acheck, hash_stats = gil.check(arr_a_truncated, show_triangle)

# gap combos 

print()
hash_gaps = {} 
nobs = len(arr_a)
for k in range(1,nobs):
    gap = arr_a[k] - arr_a[k-1]
    block = int(4*k/nobs)
    if (gap, block) in hash_gaps:
        hash_gaps[(gap, block)] += 1
    else:
        hash_gaps[(gap, block)] = 1

for gap in range(2, 16, 2):
    for block in range(4):
        key = (gap, block)
        if key in hash_gaps:
            count = hash_gaps[key]
        else:
            count = 0
        print("Gap combos:", key, count)
print("len arr_a", len(arr_a))
   

#--- 6. Impact of removing/adding one number

arr_list = [[2, 3, 5, 7, 11, 13],
            [2, 3, 5, 7, 13],
            [2, 3, 5, 7, 11, 13, 19],
            [2, 3, 5, 9, 13, 17, 25, 49, 81, 121, 169, 225, 313],
            [2, 3, 5, 9, 13, 17, 25, 49, 81, 121, 169, 225, 301, 455, 753, 1347, 2549, 2567, 2999, 3439],
            [2, 3, 5, 9, 13, 11, 15, 9, 7, 11, 21, 25, 27, 21, 13, 17]]  

print("\n--- Part 6: Removing/adding one term\n")

for arr_a in arr_list:
    flag, acheck, hash_stats = gil.check(arr_a, show_triangle)
    print("\n")


#--- 7. Continuous time series: Brownian motion

rng = np.random.default_rng(45)
np.random.seed(45)

# Generate a single Poisson random variable with lambda (lam) = 4
single_val = rng.poisson(lam=4)
arr_a = [2, 3]
arr_b = [2, 3]
a = arr_a[-1]
b = arr_a[-1]

for k in range(2000):
   sign = np.random.uniform(0, 1)
   if sign < 0.5:
       sign = -1
   else:
       sign = 1
       llambda_a = 0.8*np.log(1+k) 
       llambda_b = 0.5*np.log(1+k) 

   a = a + 2 * sign * rng.poisson(lam=llambda_a) 
   b = b + 2 * sign * rng.poisson(lam=llambda_b)
   arr_a.append(a)
   arr_b.append(b)

print("\n--- Part 7: Brownian motion\n")
show_triangle = 'Compact'
flag_a, acheck_a, hash_stats_a = gil.check(arr_a, show_triangle)
arr_stats_a = hash_stats_a['arr_amax']
flag_b, acheck_b, hash_stats_b = gil.check(arr_b, show_triangle)
arr_stats_b = hash_stats_b['arr_amax']

arr_k = np.arange(0, len(arr_a), 1)
plt.scatter(arr_k, arr_a, linewidth = 0.0, s = 1.6, c= 'red')  
plt.scatter(arr_k, arr_b, linewidth = 0.0, s = 1.6, c= 'green')

plt.gca().set_facecolor('black')
plt.show()
print(flag_a, len(arr_stats_a), flag_b, len(arr_stats_b)) 


#--- 8. Continuous time series: Smooth curve

arr_a = [2, 3, 5, 9, 13, 17]
arr_b = [2, 3, 5, 9, 13, 17]
a0 = arr_a[-1]
b0 = arr_a[-1]

for k in range(1000):
    a = a0 + 2*int(38*np.sin(k/21) + 23*np.sin(k/17)) 
    b = b0 + 2*int(28*np.sin(k/21) + 23*np.sin(k/17)) 
    arr_a.append(a)
    arr_b.append(b)

print("\n--- Part 8: Smooth curve\n")
show_triangle = 'Compact'
flag_a, acheck_a, hash_stats_a = gil.check(arr_a)  
arr_stats_a = hash_stats_a['arr_amax']
flag_b, acheck_b, hash_stats_b = gil.check(arr_b)
arr_stats_b = hash_stats_b['arr_amax']


arr_k = np.arange(0, len(arr_a), 1)
plt.scatter(arr_k, arr_a, linewidth = 0.0, s = 3.0, c= 'red') 
plt.scatter(arr_k, arr_b, linewidth = 0.0, s = 3.0, c= 'green')
plt.gca().set_facecolor('black')
plt.show()
print(flag_a, len(arr_stats_a), flag_b, len(arr_stats_b)) 


#--- 9. Build all short sequences with n <= nlevels, sigma in (0, 6, 2) 

left_start = 1
nlevels = 5

parameters = { 
               'prime': [1.0, 1.3, 1.50, 5],
               'power': [0.9, 1.1, 1.25, 5, 1.50],
             }

model = 'prime'

hash_count = {}
hash_count_valid = {}
hash_count_bounded = {}
hash_count_smooth = {}


print("\n--- Part 9: Build all short sequences starting from bottom\n")
gil.test_primes(parameters)
print()

for seed in range(0, 6, 2):
    # seed is denored as sigma in the paper

    full_hash = gil.get_all_sequences(seed, nlevels, left_start)

    for first_differences in full_hash:

        level = full_hash[first_differences]  
        key = (seed, level)
        gil.update_hash(hash_count, key, 1)
        valid, bounded, smooth, fvals, sequence = gil.admissibility(first_differences, parameters, model)
        if valid:
             gil.update_hash(hash_count_valid, key, 1)

        if valid and bounded:
            gil.update_hash(hash_count_bounded, key, 1)
            if smooth: 
                gil.update_hash(hash_count_smooth, key, 1)
                if level >= 2:
                    print("output part 9A", seed, level, sequence)

for key in hash_count:
    all = gil.get_hash(hash_count, key)
    valid = gil.get_hash(hash_count_valid, key)
    bounded = gil.get_hash(hash_count_bounded, key)
    smooth = gil.get_hash(hash_count_smooth, key)
    print("output part 9B:", key, all, valid, bounded, smooth)
