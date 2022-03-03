import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import multivariate_normal

def points_to_gaussian_heatmap(centers, height, width, scale):
    gaussians = []
    for y,x in centers:
        s = np.eye(2)*scale
        g = multivariate_normal(mean=(x,y), cov=s)
        gaussians.append(g)

    # create a grid of (x,y) coordinates at which to evaluate the kernels
    x = np.arange(0, width)
    y = np.arange(0, height)
    xx, yy = np.meshgrid(x,y)
    xxyy = np.stack([xx.ravel(), yy.ravel()]).T
    
    # evaluate kernels at grid points
    zz = sum(g.pdf(xxyy) for g in gaussians)

    img = zz.reshape((height,width))
    return img

W = 800  # width of heatmap
H = 400  # height of heatmap
SCALE = 64  # increase scale to make larger gaussians
CENTERS = [(100,100), 
           (100,300), 
           (300,100)] # center points of the gaussians

img = points_to_gaussian_heatmap(CENTERS, H, W, SCALE)

plt.imshow(img); plt.show()