import numpy as np
import gmpy2
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 'x-small'

ndigits = 10000 
ctx = gmpy2.get_context()  
ctx.precision = ndigits   
z = gmpy2.log(2)

# choose N < ndigits/2 to avoid losing precision
N = 5000  
mode = 'quantum'  # options: 'quantum' or 'chaotic'

x0 = gmpy2.mpfr(0)
x1 = gmpy2.mpfr(1)
x2 = gmpy2.mpfr(1)
arr_k = []
arr_x = []
arr_col = []
arr_log_rho = []

def set_color(k, mode):

    if mode == 'quantum':
        if k % 7 == 0:
            color='red'
        elif k % 7 == 1:
            color='blue'
        elif k % 7 == 2:
            color='green'
        elif k % 7 == 3:
            color='gold'
        elif k % 7 == 4:
            color='orange'
        elif k % 7 == 5:
            color='magenta'
        elif k % 7 == 6:
            color='gray'

    elif mode == 'chaotic':
        if k % 4 == 0:
            color = 'red'
        elif k % 4 == 1:
            color = 'gold'
        elif k % 4 == 2:
            color = 'green'
        elif k % 4 == 3:
            color = 'magenta'
        else:
            color = 'gray'
    else:
        print("Unsuported mode:", mode)
        exit()
    return(color)


#--- Main

modulus = 4
##arr_mod = np.zeros(modulus)
##arr_cnt = np.zeros(modulus)
hash_modulo = {}

for k in range(2,N):
    color = set_color(k, mode)
    if mode == 'quantum':
        x = 2*x2 - 16*x1 + 4*x0 
        delta = 0
    elif mode == 'chaotic':
        x = abs(x2 - 3*x1)
    v = gmpy2.log(abs(x))/k
    w = gmpy2.exp(v)
    rho     = float(x/x2)
    if mode == 'chaotic':   
        old_rho = float(x2/x1)
        # delta should be 0 for k < ndigits/2
        delta = rho**2 - 1 + 6/old_rho - 9/old_rho**2  
    log_rho = np.log(abs(rho))
    if k > N/10:  
        print("%6d log_rho: %8.5f %8.5f %f" % (k, log_rho, delta, rho))
        arr_k.append(k)
        arr_x.append(w)  
        arr_col.append(color)
        arr_log_rho.append(log_rho)
        residue = k % modulus
        if residue in hash_modulo:
            arr_local = hash_modulo[residue]
            arr_local.append(log_rho)
            hash_modulo[residue] = arr_local
        else:
            hash_modulo[residue] = [log_rho]
    x0 = x1
    x1 = x2
    x2 = x

plt.scatter(arr_k, arr_x, c=arr_col, s=0.02)
plt.show()
plt.scatter(arr_k, arr_log_rho, c=arr_col, s=0.02)
plt.show()
plt.plot(arr_k, arr_log_rho, linewidth=0.4)
plt.show()
print("\nRho[residue] with modulus = %2d" %(modulus))
avg = 1
for residue in range(modulus):
    arr_vals = hash_modulo[residue]
    arr_vals = np.array(arr_vals)
    arr_exp = np.exp(arr_vals)
    iqr_rho = np.quantile(arr_exp,0.75) - np.quantile(arr_exp, 0.25)
    mean_rho = np.exp(np.average(arr_vals))
    median = np.std(arr_vals)
    avg *= mean_rho
    print("residue %2d | rho = %8.5f | irq = %8.5f | median = %8.5f" 
                % (residue, mean_rho, iqr_rho, median))
avg = avg**(1/modulus)
print("Rho: %8.5f" %(avg))


#--- Plot EPDFs

np_log_rho = np.array(arr_log_rho)
import seaborn as sns
plt.figure(figsize=(8, 5))
# sns.ecdfplot(np_rho)
sns.kdeplot(data=np_log_rho, fill=True)
plt.grid(True)
plt.show()

meanlog = np.mean(np_log_rho)
medianlog = np.median(np_log_rho)
mean = np.exp(meanlog)
median = np.exp(medianlog)
print("\nRho: %8.5f [median = %8.5f | meanlog = %8.5f]" % (mean, median, meanlog))
print()


#--- Compute autocorrels for log(x/x2)

nobs = len(np_log_rho)
arr_lag = []
arr_autocorrel = []

for lag in range(200):
    mean1 = np.mean(np_log_rho[0: nobs-lag])
    mean2 = np.mean(np_log_rho[lag:nobs])
    std1 = np.std(np_log_rho[0: nobs-lag])
    std2 = np.std(np_log_rho[lag:nobs])
    dotprod = np.dot(np_log_rho[0: nobs-lag], np_log_rho[lag:nobs])
    dotprod /= (nobs-lag)
    autocorrel = (dotprod - mean1*mean2)/(std1*std2)
    arr_lag.append(lag)
    arr_autocorrel.append(autocorrel)
    print("Autcorrel lag %3d: %8.5f" %(lag, autocorrel))

plt.plot(arr_lag, arr_autocorrel, linewidth = 0.5)
plt.show()




