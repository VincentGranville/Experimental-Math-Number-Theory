import matplotlib.pyplot as plt
from matplotlib import cm # color maps
import numpy as np
import mpmath

########### ortho: plt.quiver // called direction fields
########## create maximum ascent map
############# no basin, every rain drop falling in the band ends up on the axis Re(s) = 0.5 following a smooth simple path
#### spectrum shows what happens at sigma=0.5 for a fiven t, propagates to sigma = 1 (and also below sigma = 0.5)
#######    plot between t=150 and 250
#### plot path from a rain drop
#### what other functions have this behavior? sin or sinh in complex plane??? functions with infinite hadamard product
#### can hydrology solves the most famous mathematical problem of all times??
#### plot curve at sigma = 0.5 vs sigma = 1
#### add an artificial zero or basin in the fction to see impact
#### compare with basins on earth: great basin above seal level, but not death valley!!
#### steepest descent algo .. familiar to ML practitioners
#### all minimal are the zeros
#### compare with Dirichlet eta and other fcts with infinite product : sin or (sinh z)*(sinh(sqrt 2 z)) or sinh(z sqrt(z))
#### is there always a unidirectional path?? getting to closest root? use a band b/w 2 maxima
#### sin 2x = 2 sin x cos x = ... use iteratively
#### study at sigma=3 (using derivative to find min) to find roots at sigma < 1 bc it is a lot faster// eta or zeta
### macro / micro (micro basins?) / local levels
### rescale to make spectral lines vertical, eg better with dircichlet eta than zeta
### update video to fit with za
#### fixed t, sigma varies (opposite of standard view)
### gradient on grid works with chaotic fctions
#### formula for zeta based on cosine terms
#### conformal map / synthetic function ref. to book
#### use my cosine expansion formula based on integer values to compute modulus of zeta at sigma + i t
#### title: generative AI yields novel approach to the RH, and converserly [synthetic functions]
####### https://problemsolvingwithpython.com/06-Plotting-with-Matplotlib/06.15-Quiver-and-Stream-Plots/
   ### see example using gradient
   ### https://math24.net/orthogonal-trajectories.html
   ### https://medium.com/analytics-vidhya/visualize-the-gradient-descent-of-a-cost-function-with-its-level-circles-python-d8c850731b0a
### can end up on different root depending on step size
######### grid layers with increased granularity
######### monotonic transform of za to accelerate convergence
##### https://en.wikipedia.org/wiki/Weierstrass_factorization_theorem
###       https://mathworld.wolfram.com/HadamardProduct.html
###       https://en.wikipedia.org/wiki/Riemann_zeta_function#Hadamard_product
#### moderate descent vs. steepest descent [using unequal steps x y in the grid] 
### curves with gene > 1 pick more than 1 as level not identical // start with 2 curves to create ortho

View = 'Local'       # options: 'Local' or 'Global'
Function = 'Sinh2'   # options: 'Zeta', 'Eta', 'Sinh1', 'Sinh2'
Contour = 'Lines'    # options: 'Lines', 'Surface'
Video = 'True'       # options: 'True' or 'False'
Granularity = 'Low'  # options: 'Low' or 'High'

if View == 'Local':
    min_t = 201.4    # for zeta, choose 201.0 
    max_t = 202.601  # for zeta, choose 203.301 
    min_sigma = 0.26 # for zeta, choose 0.25 
    max_sigma = 0.85 # for zeta, choose 1.60 
    if Granularity == 'High':
        incr_t = 0.0025 # slow, also requires 4x more memory
    elif Granularity == 'Low':
        incr_t = 0.005 # 4x faster than 'High', less accurate 
    incr_sigma = incr_t*(max_sigma - min_sigma)/(max_t - min_t)    
elif View == 'Global':
    min_t = 165
    max_t = 220
    min_sigma = 0.45
    max_sigma = 2
    incr_t = 0.1 
    incr_sigma = 0.01
nlevels = 180 # number of contour levels (120 for zeta)

#--- Store values of bivariate function in grid za[( , )]
#     xa[], y[a] are the x and y coordinates

xa = np.arange(min_t, max_t, incr_t)
ya = np.arange(min_sigma, max_sigma, incr_sigma)
xa, ya = np.meshgrid(xa, ya)
k_steps = 1 + int((max_t - min_t)/incr_t)
h_steps = 1 + int((max_sigma - min_sigma)/incr_sigma)
za = np.abs(0*xa + 0*ya) # set dimensions for za

k=0
for t in np.arange(min_t, max_t, incr_t):
    print("t=",t) 
    h = 0
    for sigma in np.arange(min_sigma, max_sigma, incr_sigma):    
        if Function == 'Zeta':
            z = mpmath.zeta(complex(sigma, t)) 
        elif Function == 'Eta':
            z = mpmath.altzeta(complex(sigma, t))
        elif Function == 'Sinh2':
            p = np.log(2)
            z = mpmath.cosh(complex(sigma-0.5, t*np.log(t)/4))* \
                   mpmath.cosh(complex(sigma-0.7, p*t)) 
        elif Function == 'Sinh1':
            z = mpmath.cosh(complex(sigma-0.5, t*np.log(t)/4))
        modulus=np.sqrt(z.real*z.real + z.imag*z.imag)
        za[h,k]=modulus
        h = h + 1
    k = k + 1

#--- 3D surface plot 

fig = plt.figure(figsize =(12, 8), dpi=240)
axes = plt.axes(projection ='3d')
axes.set_xlabel('t', fontsize=6)
axes.set_ylabel('sigma',fontsize=6)
axes.set_zlabel('Eta Modulus',fontsize=6)
axes.tick_params(axis='both', which='major', labelsize=6)
axes.tick_params(axis='both', which='minor', labelsize=6)
axes.set_xlim3d(min_t, max_t)
axes.set_ylim3d(min_sigma, max_sigma)
axes.set_zlim3d(0.0, float(np.max(za)))

surf = axes.plot_surface(xa, ya, za, cmap=cm.coolwarm,
                       linewidth=1, antialiased=True)

#--- Create video of 3D surface plot (rotate the plot)

if Video:

    import moviepy.video.io.ImageSequenceClip  # to produce mp4 video
    from PIL import Image  # for some basic image processing

    Nframes = 30 
    flist=[]               # list of image filenames for the video
    w, h, dpi = 4, 3, 300  # width and heigh in inches
    fps=10                 # frames per second

    for frame in range(0,Nframes): 
        image='RH4_'+str(frame)+'.png' # filename of image in current frame
        print("Creating image",image) # show progress on the screen
        axes.view_init(30, 80+ 10 * frame)
        plt.savefig(image,bbox_inches='tight')
        im = Image.open(image)
        if frame == 0:  # all images must have the same size
            width, height = im.size
            width=2*int(width/2)
            height=2*int(height/2)
            fixedSize=(width,height)
        im = im.resize(fixedSize) 
        im.save(image,"PNG")
        flist.append(image)

    # output video 
    clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(flist, fps=fps) 
    clip.write_videofile('RH4.mp4')

#--- Create contour map

fig = plt.figure(figsize =(12, 8), dpi=240)
axes = plt.axes()
CS = axes.contour(xa, ya, za, levels=nlevels, cmap=cm.coolwarm, linewidths=0.75) 
fig.colorbar(surf, ax = axes, shrink = 0.5, aspect = 5)

# Add horizontal dashed line at sigma = 0.5 (the 'critical line') 
#     https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html
plt.plot((min_t,max_t),(0.5,0.5),color='black',linestyle=(0,(15,15)),linewidth=0.2)
plt.savefig('RH4_contours.png',bbox_inches='tight')


#--- Steepest descent on contour map: sample paths to minimum (a root of |Zeta|)
#    requires: Granulary = 'High', View = 'Local'

def gradient_descent(t, sigma, showStart, showEnd, showPath, mode, n_iter, \
    learn_t, learn_sigma, type):

#   mode = Ascent or Descent, starting point is (t, sigma)
#   type = Gradient (orthogonal trajectory) or Contour
#   learn = learning rate in gradient method; n_iter = number of iterations
#   showStart, showEnd display start/end points if true, on the plot
#   in all cases, the full path is displayed in red on the plot

    x = []
    y = []
    x.append(t)
    y.append(sigma)
    h = int((sigma - min_sigma)/incr_sigma)
    k = int((t - min_t)/incr_t)
    old_z = za[h,k]
    if showStart:
        plt.plot(t, sigma, marker="o", markersize=5, markerfacecolor="yellow",color="black")
    if mode == 'Descent': 
        sign = +1
    elif mode == 'Ascent':
        sign = -1
         
    iter = 0
    while iter < n_iter:
        if type == 'Gradient':
            t = t - learn_t * sign * dx[h,k]/norm[h,k]  
            sigma = sigma - learn_sigma * sign * dy[h,k]/norm[h,k]    
        elif type == 'Contour':
            t = t - learn_t * sign * dy[h,k]/norm[h,k]  
            sigma = sigma + learn_sigma * sign * dx[h,k]/norm[h,k]    
        x.append(t)
        y.append(sigma)
        old_z = za[h, k]
        h = int((sigma - min_sigma)/incr_sigma)
        k = int((t - min_t)/incr_t)
        if h<h_steps-2 and k<k_steps-2 and h>0 and k>0: 
            z = za[h, k]
        else:
            iter = 99999999999
            showEnd = False
        iter = iter + 1

    # smooth the path x, y
    n = len(x)
    if smooth > 0:
        for i in range(1,n-2):
            x[i] = np.mean(x[max(0,i-smooth):min(n-1,i+smooth)])
            y[i] = np.mean(y[max(0,i-smooth):min(n-1,i+smooth)])

    # plot path from intitial point to convergence
    if showPath:
        plt.plot(x,y,color='red',linewidth=0.2)
    if showEnd:   
        # show where the iteration converged
        sigma=int(0.5+100*sigma)/100
        t=int(0.5+100*t)/100
        plt.plot(t, sigma, marker="o", markersize=5, markerfacecolor="palegreen",color="black")

    return(x, y)  

#--- Steepest descent: main part 

if View != 'Local':
    print('View = Local required')

dy, dx = np.gradient(za)      # matrices with same dim as za 
norm = np.sqrt(dx*dx + dy*dy) # matrix with same dim as za
angle = np.arctan2(dx,dy)     # angle of descent (unused)

smooth = 0                    # integer, to smooth tajectories (0 = no smoothing)
learn_t = incr_t              # learning parameter in gradient method  
learn_sigma = incr_sigma      # learning parameter in gradient method 

showStart = False
showEnd   = False
showPath  = False

# sample points on a contour level
#     need Ascent and Descent for full loop of contour level

level = 40 # level=1 is the one with lowest za[,]. For zeta, set level=36 
           # level must be an integer between 1 and nlevels
lines = []  
x = []
y = []
for line in CS.collections[level].get_paths():
    lines.append(line.vertices)
x = lines[0][:, 0]
y = lines[0][:, 1]
n = len(x)
showStart = False
showEnd   = False
showPath  = True
n_iter = 3000 
step = 5 # for zeta choose step=10

for i in range(0, n, step):  
    t = x[i]
    sigma = y[i]
    # need Ascent and Descent for full path up and down
    gradient_descent(t, sigma, showStart, showEnd, showPath, 'Descent', n_iter, \
       learn_t, learn_sigma, 'Gradient')
    gradient_descent(t, sigma, showStart, showEnd, showPath, 'Ascent', n_iter, \
       learn_t, learn_sigma, 'Gradient')

plt.savefig('RH4_ortho.png',bbox_inches='tight')
