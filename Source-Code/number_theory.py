from PIL import Image, ImageDraw    

n = 1000 
L = 2*n
H = int(1.1*n)

height,width = (H+1, L+1)
img1  = Image.new( mode = "RGBA", size = (L+3, H+2), color = (0, 0, 0) )
pix1  = img1.load()   
draw1 = ImageDraw.Draw(img1,"RGBA")

prod = 2**n + 1   
offset_H = 0
arr_count1 = []
arr_count0 = []
arr_count1.append(2)
arr_count0.append(n-2)


#--- 1. Main

for k in range(1, H+1): 

    prod = prod*prod  # this is (2^n + 1) at power 2^k (truncated)
    pstri = bin(prod)
    stri = pstri[0: L+2]  
    prod = int(stri, 2)
    lsr = len(stri)

    if k == n+1:
        offset_H = 1
        for l in range(L+1):
            pix1[l, k-1] = (0, 0, 0)

    stri = stri[2: len(stri)]
    if k == n:
        e_approx = stri
    rstri = stri[n:2*n]  # rightmost n digits in first 2n digits
    rcnt0 = rstri.count('0')
    rcnt1 = rstri.count('1')
    estri = stri[0:n]    # leftmost n digits
    ecnt0 = estri.count('0')
    ecnt1 = estri.count('1')
    bnum = 0

    for l in range(n):
        bnum += int(estri[l]) / 2**(l-1)
    offset_L = 0
    for l in range(min(L, len(stri))):
        if l == n+1:
            pix1[l, k-1] = (0, 0, 0)
            offset_L = 1 
        elif l == 2*n+2:
            pix1[l+1, k-1] = (0, 0, 0)
            offset_L = 2
        if stri[l] == '1':
            pix1[l+offset_L, k-1+offset_H] = (255, 0, 0) 
        elif stri[l] == '0':
            pix1[l+offset_L, k-1+offset_H] = (255, 255, 255) 

    arr_count1.append(ecnt1)
    arr_count0.append(ecnt0)

    print("%3d %3d %3d %3d %3d %f" %(k, ecnt0, ecnt1, rcnt0, rcnt1, bnum)) 

img1.save("img_1.png")
img2 = img1.crop((n-n/2, n-n/4, n, n))  # left, top, right, bottom
img2.save("img_2.png")


#--- 2. Fast computation of binary digits of e

from mpmath import mp
import numpy as np

# Set precision for n binary digits
mp.dps = int(n*np.log2(10))
e_value = mp.e  # Get e in decimal

# Convert to binary
e_binary = bin(int(e_value * (2 ** n)))[2:]

k = 0
print()
while e_approx[k] == e_binary[k]:
    k += 1
print("%d correct digits (n = %d)" %(k, n))


#--- 3. Compute nu values in Table 1

arr_rho = []
print()
for j in range(1,6):
    threshold = int(n*j/(j+1))
    if threshold % 2 == 1:
        threshold += 1 
    p1_even = arr_count1[threshold]/n
    print("Rho %1d: %7.5f" %(j, p1_even))


#--- 4. Create the plots

import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['axes.linewidth'] = 0.5
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8

xvalues = np.arange(1, H+2, 1)
arr_color = []
for k in range(H+1):
    if k%2 == 0:
        arr_color.append((1,0, 0)) 
    else:
        arr_color.append((0, 0, 1))

xmin = int(0.6*n)
xmax = int(1.1*n)    
#plt.scatter(xvalues, arr_count0, s = 0.01, c = arr_color)

plt.scatter(xvalues[xmin:xmax], arr_count1[xmin:xmax], s = 0.01, c = arr_color[xmin:xmax])
plt.axvline(x=2*n/3, color='black', linestyle='--', linewidth = 0.4, dashes=(5, 10))
plt.axvline(x=3*n/4, color='black', linestyle='--', linewidth = 0.4, dashes=(5, 10))
plt.axvline(x=4*n/5, color='black', linestyle='--', linewidth = 0.4, dashes=(5, 10))
plt.axvline(x=5*n/6, color='black', linestyle='--', linewidth = 0.4, dashes=(5, 10))
plt.axvline(x=6*n/7, color='black', linestyle='--', linewidth = 0.4, dashes=(5, 10))
plt.axvline(x=7*n/8, color='black', linestyle='--', linewidth = 0.4, dashes=(5, 10))
plt.axvline(x=8*n/9, color='black', linestyle='--', linewidth = 0.4, dashes=(5, 10))
plt.axhline(y=n/5, color='black', linestyle='--', linewidth = 0.4, dashes=(5, 10))
plt.axvline(x=n, color='red', linestyle='-', linewidth = 0.4, dashes=(5, 10))
plt.axhline(y=n/2, color='red', linestyle='--', linewidth = 0.4, dashes=(5, 10))
plt.xlim([xmin, xmax])
plt.show()
