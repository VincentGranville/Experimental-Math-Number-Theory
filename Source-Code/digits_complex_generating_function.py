import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import gmpy2

ndigits = 200
ctx = gmpy2.get_context()  
ctx.precision = ndigits

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 'x-small'
plt.figure(facecolor='black')
plt.gcf().set_facecolor("white")   # outside the plot
plt.gca().set_facecolor("black")   # plot area

xi1 = gmpy2.const_pi()
xi2 = gmpy2.exp(1)
xi3 = gmpy2.log(2)
xi4 = gmpy2.mpfr(24)/37
arr_xi = [(xi1, 'π'), (xi2, 'exp(1)'), (xi3, 'log(2)'), (xi4, '4/7')]

for m in range(len(arr_xi)): 
    value = arr_xi[m]
    xi = value[0]
    xi = gmpy2.mpz(2**(2*ndigits) * xi)
    xi_bin = bin(xi)[2:ndigits+2]
    text = value[1]
    arr_z = []
    arr_G = []
    for x in np.arange(-0.75, 0.75, 0.001):
        sum = 0
        for k in range(ndigits):
            sum += int(xi_bin[k]) * x**k 
        arr_z.append(x)
        arr_G.append(sum)
    plt.plot(arr_z, arr_G, linewidth=1.4, label = f"ξ = {text}")
plt.grid(True, color = (0.3, 0.3, 0.3))
plt.legend(loc='upper left', shadow=True)
plt.show()  
 

#--- G(z) orbit on the unit circle

plt.gcf().set_facecolor("white")   # outside the plot
plt.gca().set_facecolor("black")   # plot area

xi = gmpy2.mpz(2**(2*ndigits) * xi1)  
xi_bin = bin(xi)[2:ndigits+2]

## arr_z = []
arr_G_re = []
arr_G_im = []
arr_colors = []
npoints = 200000  # points generated on the circle
radius = 0.9999995
for t in range(npoints):
    tc = t/npoints
    color = [tc, 4*tc*(1-tc), 0.25+0.75*(1-tc)]
    z_re = radius * np.cos(2*np.pi*tc)
    z_im = radius * np.sin(2*np.pi*tc)
    if t % 1000 == 0:
        print("...",t)
    sum = 0
    for k in range(ndigits):
        sum += int(xi_bin[k]) * complex(z_re, z_im)**k 
    sum /= -np.log(1 - radius)
    if abs(sum.real) < 0.5 and abs(sum.imag) < 0.5:
        arr_G_re.append(sum.real)
        arr_G_im.append(sum.imag)
        arr_colors.append(color)

plt.scatter(arr_G_re, arr_G_im, s=0.4, color = arr_colors, edgecolors = 'none')
plt.show()
