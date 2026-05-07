# prng_test_vg.py 

import numpy as np
import lzma 


def test_frequencies(Lmax, rnd_bits):   
    
    M = len(rnd_bits) 
    print("\nString", end = " ")
    for test in range(M):
        print("  Test %1d" % test, end = " ")
    print()
    for L in range(1,Lmax+1):
        for k in range(2**L):
            arr_stats = []
            kbin = bin(k)[2:]
            kbin = '0' * (L-len(kbin)) + kbin
            for m in range(M):
                bits = rnd_bits[m]
                indices = [i for i in range(len(bits)) if bits.startswith(kbin, i)]
                arr_stats.append(len(indices))
            print("%6s" % (kbin), end = " ") 
            for test in range(len(arr_stats)):
                print("%8d" % (arr_stats[test]), end=" ")
            print()
    return()


def test_autocorrel(maxLag, rnd_bits):  

    M = len(rnd_bits)  
    length = min(len(rnd_bits[0])-maxLag, 5000000)
    arr_max_correl = []

    for m in range(M):
        bits = rnd_bits[m]
        auto_correl = []
        for lag in range(1, maxLag):
            arr1 = np.array([int(b) for b in bits[0:length]])
            arr2 = np.array([int(b) for b in bits[lag:lag+length]])
            correl = np.corrcoef(arr1, arr2)[0, 1]
            auto_correl.append(correl) 
        auto_correl = np.array(auto_correl)
        max_ac = np.max(np.abs(auto_correl))
        arr_max_correl.append(max_ac)

    print("\nMax abs autocorrels up to lag %2d:" %(maxLag-1))
    for test in range(len(arr_max_correl)):
        print("Test %2d: %8.6f" % (test, arr_max_correl[test]))
    return()


def update_hash(hash, key, count):
    if key in hash:
        hash[key] += count
    else:
        hash[key] = count
    return(hash)


def test_collisions(L, rnd_bits): 

    M = len(rnd_bits)  
    length = len(rnd_bits[0]) - L
    arr_collisions = []
 
    for m in range(M):
        bits = rnd_bits[m]
        hash = {}
        for idx in range(0, length, 1): 
            update_hash(hash, bits[idx:idx+L], 1)
        collisions = 0
        max_hits = 0
        for string in hash:
            cnt = hash[string]
            if cnt > 1:
                collisions += 1
            if cnt > max_hits:
                max_hits = cnt
        arr_collisions.append([len(hash), collisions, max_hits])

    print("\n%2d-bits   collisions max_hits" %(L))
    for test in range(len(arr_collisions)):
        collisions = arr_collisions[test]
        print("Test %2d:   %8d %5d" % (test, collisions[1], collisions[2]))
    return()


def test_spectral(rnd_bits): 
    
    from math import erfc 
    arr_spectral = [] 
    M = len(rnd_bits)  

    for m in range(M):

        bits = rnd_bits[m]
        x = np.array([1 if b == '1' else -1 for b in bits], dtype=float)
        n = len(x)
        s = np.fft.fft(x)
        m = np.abs(s[:n // 2])
        tau = np.sqrt(np.log(1.0 / 0.05) * n)  # threshold
        n0 = 0.95 * n / 2.0    # expected peaks below threshold
        n1 = np.sum(m < tau)   # peaks below threshold
        d = (n1 - n0) / np.sqrt(n * 0.95 * 0.05 / 4.0)  # test statistic
        p_value = erfc(abs(d) / np.sqrt(2.0))
        result = 'FAIL'
        if p_value >= 0.01:
            result = 'PASS'
        arr_spectral.append([p_value, result])
    
    print("\nspectral    p_value")
    for test in range(M):
        spectral = arr_spectral[test]
        print("Test %2d:   %8.5f   %s" % (test, spectral[0], spectral[1]))
    return()


def test_compress(rnd_bits): 

    M = len(rnd_bits)  
    arr_compress = []
 
    for m in range(M):
        bit_string = rnd_bits[m]

        padded_bit_string = bit_string.ljust((len(bit_string) + 7) // 8 * 8, '0')
        byte_data = int(padded_bit_string, 2).to_bytes(len(padded_bit_string) // 8, byteorder='big')
        compressed_data = lzma.compress(byte_data, preset=9)
        arr_compress.append(len(compressed_data))

    print("\nLZMA       size   ratio")
    for test in range(len(arr_compress)):
        size = arr_compress[test]
        ratio = size / (len(bit_string)/8) 
        print("Test %2d: %6d  %8.5f" % (test, size, ratio))
    return()


def test_runs(bit, rnd_bits):  

    import matplotlib.pyplot as plt
    import matplotlib as mpl
    M = len(rnd_bits)  

    mpl.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.labelsize'] = 8
    plt.rcParams['ytick.labelsize'] = 8
    plt.rcParams['legend.fontsize'] = 'x-small'

    length = len(rnd_bits[0])
    qualitative_colors = plt.get_cmap('tab10').colors
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))
    ax1.set_facecolor('black')
    ax2.set_facecolor('black')

    for m in range(M):

        if m == 0:
            color = 'blue'
        else:
            color = qualitative_colors[m+1]

        arr_run = []
        arr_run_freq = []
        arr_run_first = []

        bit_string = rnd_bits[m]
        hash_runs = {}
        hash_records = {}
        run = 1
        for t in range(1, length):
            Flag = True
            if bit_string[t] == bit_string[t-1] == str(bit):
                run += 1
            else:
                Flag = False
                if run not in hash_runs:
                    hash_records[run] = t
                update_hash(hash_runs, run, 1)
                run = 1
                
        for run in hash_runs:
            if run > 1:
                arr_run.append(run)
                arr_run_freq.append(np.log(hash_runs[run]))
                arr_run_first.append(np.log(hash_records[run]))
 
        ax1.plot(arr_run, arr_run_freq, color=color, linewidth = 0.8, label=f"Test {m}")
        ax2.plot(arr_run, arr_run_first, color=color, linewidth = 0.8, label=f"Test {m}")

    ax1.legend(loc='upper right')
    ax2.legend(loc='upper left')
    plt.show()
    return()


def test_next_bit(rnd_bits):  

    from sklearn.neural_network import MLPClassifier
    from scipy.stats import binom

    M = len(rnd_bits)  
    n = 1000               # samples to extract from string
    window = 10            # predict next bit based on previous bits in window
    split = int(0.80 * n)  # 80% for training, 20% for testing
    nobs = n - split
    arr_accuracy = []
    length = len(rnd_bits[0])
    np.random.seed(56)     # PCG64 has issues if this seed is 56 or 58
    offsets = np.random.randint(0, length-window, size=n)

    for m in range(M):

        bits = rnd_bits[m]
        ibits = np.array([1 if b == '1' else 0 for b in bits], dtype=int)
        # predict bit at position start+windown based on past window bits
        X, y = [], []
        for k in range(n):
            start = offsets[k]
            X.append(ibits[start : start + window])
            y.append(ibits[start + window])
        X_train, X_test = X[:split], X[split:]
        y_train, y_test = y[:split], y[split:]

        clf = MLPClassifier(hidden_layer_sizes=(32, 16), max_iter=2000, random_state=42)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        arr_errors = y_pred - y_test     # 0 if match; are errors randomly spread or not?
        accuracy = clf.score(X_test, y_test)
        arr_accuracy.append([accuracy, arr_errors])

    print("\nnext-bit accuracy p_value (0.50 accuracy means cannot predict next bit)")
    for test in range(len(arr_accuracy)):
        arr_output = arr_accuracy[test]
        accuracy = arr_output[0]
        lower = min(accuracy, 1-accuracy) * nobs
        upper = max(accuracy, 1-accuracy) * nobs
        p_valueL = binom.cdf(lower, nobs, 0.50)
        p_valueR = 1- binom.cdf(upper, nobs, 0.50)
        p_value = p_valueL + p_valueR
        print("Test %2d %8.5f %8.5f" % (test, accuracy, p_value))
    return()
