import numpy
import pyDOE
import scipy.stats
from scipy.stats.distributions import norm
import scipy.optimize as optimize
from numpy.linalg import eig, inv
from pyquaternion import Quaternion
from os.path import isfile
import glob
import multiprocessing as mp
import sys
import matplotlib.pyplot as plt

def load_normalised_voxel_growth_space(fname,resolution):
	# Function reads growth_space file, shifts x,y coordinates
    # to zero mean, sets min z coordinate to 0.
    gs = numpy.genfromtxt(fname,delimiter=',')
    gs = gs*resolution/100.
    gs[:,0] = gs[:,0] - numpy.mean(gs[:,0])
    gs[:,1] = gs[:,1] - numpy.mean(gs[:,1])
    gs[:,2] = gs[:,2] - numpy.min(gs[:,2])
    return gs

def fit_ellipse(x,y):
# Fit_ellipse function takes in 2D points in x,y (i.e. 1D arrays)
# and returns orientation of the ellipse (etheta), major and minor
# axes of the ellipse (exr1, exr2)
    # Find ellipse given points in x,y
    x = x[:,numpy.newaxis]
    y = y[:,numpy.newaxis]
    D =  numpy.hstack((x*x, x*y, y*y, x, y, numpy.ones_like(x)))
    S = numpy.dot(D.T,D)
    C = numpy.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(numpy.dot(inv(S), C))
    n = numpy.argmax(numpy.abs(E))
    a = V[:,n] # a is ellipse
    # Ellipse centre
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    ex,ey = numpy.array([x0,y0]) # Ex,Ey are coords of ellipse centre
    etheta = 0.5*numpy.arctan(2*b/(a-c)) # Ellipse angle of rotation
    # Ellipse major/minor axes
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*numpy.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*numpy.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=numpy.sqrt(numpy.abs(up/down1))
    res2=numpy.sqrt(numpy.abs(up/down2))
    exr1,exr2 = numpy.array([res1, res2])
    return etheta,exr1,exr2

def calc_growth_space_vx(pts):
    # Bounding box vals
    xmin,xmax = numpy.min(pts[:,0]),numpy.max(pts[:,0])
    ymin,ymax = numpy.min(pts[:,1]),numpy.max(pts[:,1])
    zmin,zmax = numpy.min(pts[:,2]),numpy.max(pts[:,2])
    bbox = numpy.asarray([xmin,xmax,ymin,xmax,zmin,zmax])
    # Estimation of trunk height
    r = []
    for i in range(100):
        lh,uh = (zmax/100.)*i, (zmax/100.)*i+1
        p_a_h = numpy.logical_and(pts[:,2]>lh,pts[:,2]<uh)
        sum_p_a_h = numpy.sum(numpy.logical_and(pts[:,2]>lh,pts[:,2]<uh))
        if sum_p_a_h>1:
            temp_thet,temp_r1,temp_r2 = fit_ellipse(pts[p_a_h,0],pts[p_a_h,1])
            if ((temp_r1+temp_r2)/2.<30.):                
                r.append([lh, (temp_r1+temp_r2)/2.])
    
    bottom_rad = numpy.mean(r[0:8])
    for j in range(len(r)):
        if r[j][1]>3*bottom_rad:
            trunk_height = r[j][0]
            crown_height = zmax-trunk_height;
            break
    return bbox, trunk_height, crown_height

def main():
    gs_dir = "../../growth-space/voxel_size_tests/"
    gs_filenames = ["Tree1_25cm_5pts.csv", "Tree1_50cm_5pts.csv",\
        "Tree1_75cm_5pts.csv", "Tree1_100cm_5pts.csv"]
    resolution = [10.,20.,25.,50.,75.,100.]
    print "Resolution\tTrunk height\tCrown height\tTree height"
    print "--------------------------------------------------------------------"
    for i in resolution:
        gs_filename = "Tree1_" + str(int(i)) + "cm_5pts.csv"
        tpts = load_normalised_voxel_growth_space(gs_dir+gs_filename,i)
        targ_bb, targ_th,targ_ch = calc_growth_space_vx(tpts)
        print "%s\t\t%2.2f\t\t%2.2f\t\t%2.2f" % (i,targ_th,targ_ch, targ_bb[5])

    
    
    
    
    
if __name__ == "__main__":
    main()