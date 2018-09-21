import numpy
from numpy.linalg import eig, inv


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

def cgs_ls(pts):
    # Bounding box vals
    xmin,xmax = numpy.min(pts[:,0]),numpy.max(pts[:,0])
    ymin,ymax = numpy.min(pts[:,1]),numpy.max(pts[:,1])
    zmin,zmax = numpy.min(pts[:,2]),numpy.max(pts[:,2])
    bbox = numpy.asarray([xmin,xmax,ymin,xmax,zmin,zmax])
    # Estimation of trunk height
    th1 = numpy.max(pts[numpy.logical_and(pts[:,0]==0,pts[:,1]==0),2])
    th2 = numpy.min(pts[numpy.logical_and(pts[:,0]!=0,pts[:,1]!=0),2])
    trunk_height = numpy.min([th1,th2])
    # Crown height
    crown_height = zmax-trunk_height
    # Geometric mean of crown
    mu_x = numpy.mean(pts[pts[:,2]>trunk_height,0])
    mu_y = numpy.mean(pts[pts[:,2]>trunk_height,1])
    mu_z = numpy.mean(pts[pts[:,2]>trunk_height,2])
    mu_c = numpy.array([mu_x,mu_y,mu_z])
    # Calculate radii of thirds of crown
    # Lower
    l_et1,l_e1,l_e2 = \
        fit_ellipse(pts[numpy.logical_and(pts[:,2]>trunk_height,pts[:,2]<mu_z),0],\
                    pts[numpy.logical_and(pts[:,2]>trunk_height,pts[:,2]<mu_z),1])
    el = numpy.array([l_et1,l_e1,l_e2])
    # Upper
    u_et1,u_e1,u_e2 = \
        fit_ellipse(pts[pts[:,2]>mu_z,0],pts[pts[:,2]>mu_z,1])
    eu = numpy.array([u_et1,u_e1,u_e2])
    return bbox,trunk_height,crown_height,mu_c,el,eu

def cgs_vx(pts):
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
            try:
                temp_thet,temp_r1,temp_r2 = fit_ellipse(pts[p_a_h,0],pts[p_a_h,1])
                if ((temp_r1+temp_r2)/2.<30.):
                    r.append([lh, (temp_r1+temp_r2)/2.])
            except:
                continue
    bottom_rad = numpy.mean(r[0:2])
    for j in range(len(r)):
        if r[j][1]>1.5*bottom_rad:
            trunk_height = r[j][0]
            crown_height = zmax-trunk_height;
            break
    # Geometric mean of crown
    mu_x = numpy.mean(pts[pts[:,2]>trunk_height,0])
    mu_y = numpy.mean(pts[pts[:,2]>trunk_height,1])
    mu_z = numpy.mean(pts[pts[:,2]>trunk_height,2])
    mu_c = numpy.array([mu_x,mu_y,mu_z])
    # Calculate radii of thirds of crown
    # Lower
    l_et1,l_e1,l_e2 = \
        fit_ellipse(pts[numpy.logical_and(pts[:,2]>trunk_height,pts[:,2]<mu_z),0],\
                    pts[numpy.logical_and(pts[:,2]>trunk_height,pts[:,2]<mu_z),1])
    el = numpy.array([l_et1,l_e1,l_e2])
    # Upper
    u_et1,u_e1,u_e2 = \
        fit_ellipse(pts[pts[:,2]>mu_z,0],pts[pts[:,2]>mu_z,1])
    eu = numpy.array([u_et1,u_e1,u_e2])
    return bbox, trunk_height, crown_height, mu_c, el, eu#bbox,trunk_height,crown_height,mu_c,el,eu
