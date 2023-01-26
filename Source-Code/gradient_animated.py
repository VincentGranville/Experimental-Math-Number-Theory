# replace the main part in gradient.py by the following

#--- Steepest descent: main part 
# import moviepy.video.io.ImageSequenceClip  

dy, dx = np.gradient(za)      # matrices with same dim as za 
smooth = 0                    # integer, to smooth tajectories (0 = no smoothing)
learn_t = incr_t              # learning parameter in gradient method  
learn_sigma = incr_sigma      # learning parameter in gradient method 

n = 100
np.random.seed(101)
xx = np.random.uniform(min_t,max_t,n)
yy = np.random.uniform(min_sigma,max_sigma,n)
showEnd   = False
showPath  = True
n_iter = 1 
flist = []
fps = 4

for frame in range(200):
    
    print("frame",frame)
    image='RH4_ortho'+str(frame)+'.png'
    if frame == 0:
        showStart = True
    else:
        showStart = False

    for i in range(0, n):  
        t = xx[i]
        sigma = yy[i]
        (x, y) = gradient_descent(t, sigma, showStart, showEnd, showPath, 'Descent', n_iter, \
           learn_t, learn_sigma, 'Gradient')
        m = len(x)
        xx[i] = x[m-1]  # new t attached to starting point i
        yy[i] = y[m-1]  # new sigma attached to starting point i
    plt.savefig(image,bbox_inches='tight')
    flist.append(image)

# output video 
clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(flist, fps=fps) 
clip.write_videofile('RH4_ortho.mp4'
