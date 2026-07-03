import numpy as np
import copy

#--- Utils

def make_readable(data, type = "int"):

    readable_data = []
    for elt in data:
        if type == "int":
            readable_data.append(int(elt))
        elif type == "float":
            readable_data.append(float(elt))
    return(readable_data)


def update_hash(hash, key, count):
    if key in hash:
        hash[key] += count
    else:
        hash[key] = count
    return(hash)


def get_hash(hash, key):
    if key in hash:
        count = hash[key]
    else:
        count = 0
    return(count)


def arr_plots(arr_k, arr_min, arr_max, arr_a, resets_list = ()):

    import matplotlib.pyplot as plt
    import matplotlib as mpl

    mpl.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.labelsize'] = 8
    plt.rcParams['ytick.labelsize'] = 8
    plt.rcParams['legend.fontsize'] = 'x-small'

    arr_k = np.array(arr_k)
    arr_min = np.array(arr_min)
    arr_max = np.array(arr_max)
    arr_a = np.array(arr_a)

    from matplotlib.ticker import FormatStrFormatter  
    from matplotlib.ticker import MaxNLocator

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(0.8*9, 0.8*3))

    ax1.yaxis.get_offset_text().set_fontsize(6)
    ax1.plot(arr_k, arr_min, lw=0.4, color='lightgreen')
    ax1.plot(arr_k, arr_max, lw=0.4, color='orange')
    ax1.plot(arr_k, arr_a, color='white', lw = 0.4)
    for k in resets_list:
        ax1.axvline(x=k, color='red', lw=0.4)
    ax1.ticklabel_format(style='sci', scilimits=(0, 0), axis='y', useMathText=True)
    ax1.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
    ax1.set_facecolor('black')
    ax1.tick_params(labelsize=6, pad=2, length=2)  

    ax2.plot(arr_k, arr_min/arr_a, lw=0.4, color='lightgreen') 
    ax2.plot(arr_k, arr_max/arr_a, lw=0.4, color='orange') 
    ax2.axhline(y=1.0, color='white', lw=0.4) 
    ax2.ticklabel_format(style='sci', scilimits=(0, 0), axis='y', useMathText=True)
    ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax2.yaxis.set_major_locator(MaxNLocator(nbins=5, prune='lower'))
    ax2.set_facecolor('black')
    ax2.tick_params(labelsize=6, pad=2, length=2) 
    ax2.set_ylim(0.2, 1.8)

    ax3.yaxis.get_offset_text().set_fontsize(6)
    ax3.plot(arr_k, arr_min-arr_a, lw=0.4, color='lightgreen')
    ax3.plot(arr_k, arr_max-arr_a, lw=0.4, color='orange')
    ax3.axhline(y=0.0, color='white', lw=0.4) 
    ax3.ticklabel_format(style='sci', scilimits=(0, 0), axis='y', useMathText=True)
    ax3.yaxis.set_major_locator(MaxNLocator(nbins=6, prune='lower'))
    ax3.set_facecolor('black')
    ax3.tick_params(labelsize=6, pad=2, length=2)

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.15) 
    plt.show()
    return()


#--- banned prime gaps series and simulation

def banned_gaps_loop(gaps, start):

    val = start 
    seq = (val,)
    for k in range(0, len(gaps)):
        val += gaps[k]
        seq = (*seq, val) 

    banned = False

    for modulo in (3, 5, 7, 11, 13): 
        hash_residues = {}
        for val in seq:
            residue = val % modulo
            update_hash(hash_residues, residue, 1) 
        
        if len(hash_residues) == modulo:
            banned = True
            break

    return(banned)


def banned_gaps(gaps):

    banned = True
    for start in (2, 4, 6, 8, 10, 12):
        banned_start = banned_gaps_loop(gaps, start)
        if not banned_start:
            banned = False
            break
    return(banned)


def banned_found_in_seq(first_diff_seq):

    banned_found = False
    banned_pattern = ()

    for size in (2, 3, 4, 5, 6, 7):
        for k in range(len(first_diff_seq)-size):
            gaps = ()
            for j in range(size):
                gaps = (*gaps, first_diff_seq[k+j]) 
            banned = banned_gaps(gaps)
            if banned:
                banned_found = True
                banned_pattern = gaps
                break
        if banned_found:
            break
    return(banned_found, banned_pattern)


def generate_random_sequence(n, arr_init, llambda, parameters, model, rng, rejection_sampling = True): 

    # the generated sequence is guaranteed to be valid (no duplicates)

    q = arr_init[-1]  
    arr_a = np.copy(arr_init)
    sum = np.sum(arr_a)
    offset = len(arr_init)

    for k in range(n):  
      
        ccontinue = True

        while ccontinue: 
            # until in corridor
            
            u = 2 + 2*rng.poisson(lam=llambda)  # try increasing llambda
            new_q = q + u
            current_n = offset + k + 1

            admissible = admissibility_local(q, new_q, sum, current_n, parameters, model)
            if admissible or not rejection_sampling:
                ccontinue = False

        q += u
        arr_a = np.append(arr_a, q)

    return(arr_a)


#--- get canonical form of sequence arr_a

def canonical(s1, n_levels, fail_level = {}, valid = False, threshold = 1.0):

    # s1 is in the input sequence
    # s2 is the canonical form of s1 if nor failures added during construction
    # fail_level = {} for no failure
    # fail_level = {23:3, 41:5} to build fails at level 23, 41 resp. with 3, 5 instead of 1 (1 = no fail) 

    arr_a = np.copy(s1)
    np.random.seed(66)
    n = len(arr_a)
    arr_p =[-1]

    for level in range(1, n_levels+1, 1):

        arr1 = np.copy(arr_a)
        arr_a = []
        cnt = 0
        for k in range(len(arr1) - 1):
            a = abs(arr1[k+1] - arr1[k])
            arr_a.append(a)
            if arr_a[k] < arr1[k]:
                cnt += 1
        arr_a = make_readable(arr_a)
        p = cnt / len(arr_a)
        arr_p.append(p)
        if level == 1:
            s1_max_gap = np.max(arr_a)

    s1_bottom = np.copy(arr_a)

    for level in range(n_levels-1, -1, -1):

        # "reduced" primes is last arr_a produced (after adding 1)

        arr1 = np.copy(arr_a)
        if level in fail_level:
            fail_value = fail_level[level]
            arr_a = [fail_value]
        else:
            arr_a = [1]   

        for k in range(1, len(arr1)+1):

            a1 = arr_a[k-1] + arr1[k-1]
            a2 = arr_a[k-1] - arr1[k-1]   
            if a2 < 0:
                a = a1
            else:
                u = np.random.uniform(0,1)
                # random choice
                if u < threshold:  # try u < arr_p[level]
                    a = a2
                else:
                    a = a1
            if valid:
                if level == 0:
                    # make sure next prime is not less than previous prime
                    a = a1
                elif level == 1:
                    # make sure next prime is bigger than previous one
                    if a2 == 0:
                        a = a1
            arr_a.append(a)

        arr_a = make_readable(arr_a)
        if level == 1:
            s2_gaps = np.copy(arr_a)

    for k in range(len(arr_a)):
        arr_a[k] += 1
    s2 = np.copy(arr_a)
    s1_proba = make_readable(arr_p, type = "float")

    return(s2, s2_gaps, s1_bottom, s1_max_gap, s1_proba)


#--- find bounds for next prime

def augment(right_diagonal, a):

    new_diagonal = [a,]
    val = a
    success = False
    for j in range(len(right_diagonal)):  
        next = abs(val - right_diagonal[j])
        new_diagonal.append(next)
        val = next
    if val == 1:
        success = True
    return(success, new_diagonal)


def max_gap_allowed(right_diagonal, a):

    lower = a
    upper = 2*a + 1
    found = False
    trials = 0

    while not found:
        trials += 1
        b = (upper + lower) // 2 
        if b % 2 == 0:
            b += 1
        success, new_diagonal = augment(right_diagonal, b)
        if success:
            lower = b
        else:
            upper = b 
        if upper - lower == 2 or trials > 28:
            found = True

    return(b, trials)


def min_gap_allowed(right_diagonal, a):

    lower = 1
    upper = a
    found = False
    trials = 0

    while not found:
        trials += 1
        b = (upper + lower) // 2 
        if b % 2 == 0:
            b += 1
        success, new_diagonal = augment(right_diagonal, b)
        if success:
            upper = b
        else:
            lower = b 
        if upper - lower == 2 or trials > 280:
            found = True
        
    return(b, trials)


def test_conjecture(old_a, max_a, right_diagonal):

    satisfied = True
    flag = 'no issue'

    # at 'a = max + 2', we must fail, it's above the max limit
    success1, test_diagonal = augment(right_diagonal, max_a + 2)
    success1 = bool(not success1)

    if success1:
        # check values below max_a; they must all succeed
        for a in range(old_a, max_a, 2):  
            success2, test_diagonal = augment(right_diagonal, a)
            if not success2:
                flag = (old_a, a, max_a)
                satisfied = False
                break
    else: 
        flag = 'a_max + 2 succeeds, not supposed to'
        satisfied = False
    return(satisfied, flag)


def collect_statistics(right_diagonal):

    length = len(right_diagonal)
    idx = -1
    idx2 = -1
    end = -1
    reversal_idx = - 1
    reversal_value = -1
    reversal = False
    reset = False
    hash_output = {}

    for kx in range(length-1, 0, -1):
        if right_diagonal[kx] > 2:
            idx = kx
            end = right_diagonal[kx]
            break
    if idx == -1:
        reset = True 

    dmax2 = np.max(right_diagonal[1:])
    for kx in range(length-1, 0, -1):
        if right_diagonal[kx] == dmax2:
            dmax2_idx = kx
            break

    for kx in range(len(right_diagonal[1:])):
        d1 = right_diagonal[kx]
        d2 = right_diagonal[kx+1]
        if d2 >= d1 and d2 > 2:
            reversal_idx = kx + 1
            reversal_value = d2
            reversal = True
            break

    # dmax is maximum value after the first log(length)
    start = int(np.log(length))
    dmax = np.max(right_diagonal[start:]) 

    counts = ()
    for kx in (0, 2):
        counts = (*counts, right_diagonal.count(kx))
    
    arr_cycle02 =right_diagonal[idx+1:-2]
    nu_2 = arr_cycle02.count(2)

    hash_output['reset'] = reset
    hash_output['counts'] = counts
    hash_output['max_value'] = dmax2
    hash_output['max_index'] = dmax2_idx
    hash_output['before_02_value'] = end
    hash_output['before_02_index'] = idx
    hash_output['max_tail_value'] = dmax
    hash_output['reversal_idx'] = reversal_idx
    hash_output['reversal_value'] = reversal_value
    hash_output['reversal'] = reversal
    hash_output['nu_2'] = nu_2

    return(hash_output)


#--- Main function

def check(arr_a, show_triangle = 'None', mode = 'Fast', delta = 'standard'):

    hash_stats = {}
    right_diagonal = []
    left_diagonal = []
    arr_amax = []
    right_diagonal.append(arr_a[-1])
    arr1 = np.copy(arr_a)
    arr1 = np.array(arr1).astype(int)
    if show_triangle in ('Full',): 
        print("Sequence:", arr_a) 
    imax = len(arr_a) - 1
    flag = 'Success'
    level = 0
    failing_level = -1

    for k in range(imax, 0, -1):

        level += 1
        arr2 = np.zeros(k)
        for j in range(k):
            if delta == 'standard':
                arr2[j] = abs(arr1[j+1] - arr1[j])
            elif delta == 'light':
                arr2[j] = max(abs(arr1[j+1] - arr1[j]), 1)
        if k == imax:
            first_diff = np.copy(arr2) 
        right_diagonal.append(arr2[-1])
        left_diagonal.append(arr2[0])

        amax = np.max(arr2)  
        first_amax_idx = 1 + np.argmax(arr2 == amax) 
        first_amax_val = int(arr1[first_amax_idx - 1])
        amax = int(amax)
        arr1 = np.copy(arr2)
        arr1 = np.array(arr1).astype(int)

        count = np.count_nonzero(arr1 == amax)
        arr_amax.append([level, amax, count]) 
        if arr1[0] != 1:
            flag = 'Fail'
            failing_level = level 
            if not mode == 'Full':
                break
        if show_triangle == 'Compact':
            print("Level",level,"|amax=",amax,"|idx=",first_amax_idx,"|val=",first_amax_val,"|count=",count) 
        elif show_triangle == 'Full':
            print("Level",level,"|amax=",amax,"|idx=",first_amax_idx,"|val=",first_amax_val,arr1) 
        if amax == 2:
            flag = 'Success'
            if mode == 'Fast':
                break

    right_diagonal = np.array(right_diagonal).astype(int)
    left_diagonal = np.array(left_diagonal).astype(int)
    first_diff = np.array(first_diff).astype(int)

    hash_stats['arr_amax'] = arr_amax
    hash_stats['first_diff'] = first_diff
    hash_stats['right_diagonal'] = right_diagonal
    hash_stats['left_diagonal'] = left_diagonal
    hash_stats['failing_level'] = failing_level

    return(flag, first_amax_idx, hash_stats)


def admissibility_local(q, proposed_next_q, sum, n, parameters, model):

    # testing a new value; assumes it is admissible prior to that

    constants = parameters[model]
    lower_constant = constants[0]  
    upper_constant = constants[1]
    smooth_constant = constants[2]

    old_sum = sum
    sum += proposed_next_q

    if model == 'power':
        exponent = constants[4]
        growth = n**exponent 
    elif model == 'prime':
        growth = n*np.log(n)
    lower = lower_constant * growth
    upper = upper_constant * growth 

    pass_test1 = bool(lower < proposed_next_q < upper) 
    pass_test2 = bool(proposed_next_q < smooth_constant * q) 

    if pass_test1 and pass_test2:
        admissible = True
    else:
        admissible = False
    return(admissible)


def admissibility(first_differences, parameters, model): 

    # check if sequence is in corridor; test the whole sequence
    # n = index of first failure
    # fvals = (n lower bound, observed value, upper bound) if failing bound test
    # fvals = (n, observed value, max allowed given past value) if failing smooth test

    valid = True
    bounded = True
    smooth  = True
    fvals = ()

    if 0 in first_differences:

        valid = False
        sequence = ()

    else:

        constants = parameters[model]
        lower_constant = constants[0]  
        upper_constant = constants[1]
        smooth_constant = constants[2]
        offset = constants[3]

        sum = 2
        sequence = (sum,)
        n = 1

        for value in first_differences:

            n += 1
            old_sum = sum
            sum += value
            if model == 'power':
                exponent = constants[4]
                growth = n**exponent 
            elif model == 'prime':
                growth = n*np.log(n)
            lower = lower_constant * growth
            upper = upper_constant * growth 

            pass_test1 = bool(lower < sum < upper)
            if not pass_test1 and n > offset: 
                bounded = False
                if len(fvals) == 0:
                    fvals = (n, lower, sum, upper)
            pass_test2 = bool(sum < smooth_constant * old_sum)
            if not pass_test2 and n > offset:
                smooth = False
                if len(fvals) == 0:
                    fvals = (n, lower, smooth_constant * old_sum)
            sequence = (*sequence, sum)

    return(valid, bounded, smooth, fvals, sequence)


def test_primes(parameters):

    # check if prime number sequence (first 1000 primes) is in corridor
    from primePy import primes
    prime_list = primes.first(1000)
    # prime_list = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31) 
 
    model = 'prime'
    first_differences = ()
    for k in range(len(prime_list) - 1):
        delta = prime_list[k+1] - prime_list[k] 
        first_differences = (*first_differences, delta)
    valid, bounded, smooth, fvals, sequence = admissibility(first_differences, parameters, model)
    print("Testing prime sequence: bounded = %s | smooth = %s " %(bounded, smooth))
    return()


def integrate(x, left_start = 1):

    # find all sequences y one level above x such that diff(y) = x

    y = (left_start,)
    list2 = (y,)

    for k in range(1, 1+len(x)):

        list3 = ()  # used to be {}

        for y in list2:  

            # we must have if len(y) == k:
    
            y_k_plus = y[k-1] + x[k-1]
            new_y_plus = (*y, y_k_plus)
            list3 = (*list3, new_y_plus)

            if x[k-1] !=0: 
                y_k_minus = y[k-1] - x[k-1]
                if y_k_minus >= 0:
                    new_y_minus = (*y, y_k_minus)
                    list3 = (*list3, new_y_minus)

        list2 = copy.deepcopy(list3)

    return(list2)


def get_all_sequences(seed, nlevels, left_start = 1):

    hash_arr = { (left_start, seed):1 }
    full_hash = {}

    for level in range(nlevels):

        hash3 = {}
        count = len(hash_arr)
        sub_count = 0

        for arr1 in hash_arr:

            if sub_count % 100 == 0 and sub_count > 0:
                print("GBL Progress:",seed, level, sub_count,"/", count) 
            sub_count += 1
            local_list = integrate(arr1, left_start)

            # Add the collected local_lists to list2 
            for arr2 in local_list:
                hash3[arr2] = level + 1 

        hash_arr = hash3.copy()
        full_hash.update(hash_arr)  

    return(full_hash)

