import numpy
import pyDOE
import scipy.stats
import scipy.optimize as optimize
from numpy.linalg import eig, inv
from pyquaternion import Quaternion
from os.path import isfile


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

def read_growth_space(pts):
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

def load_target_points():
    tpts = numpy.genfromtxt(target_filename,delimiter=' ',usecols=(1,2,3))
    return tpts


def estimate_error(params):
# Estimate_error function takes in:
# params: parameters used to create tree from L-System
# fname: filename of obj file (comparison data)
    tpts = load_target_points()
    targ_bbox,targ_th,targ_ch,targ_mu,targ_el,targ_eu = read_growth_space(tpts)
	# Generate L-system points for params
    ls_pts = generate_lsystem_tree_points(params)
    ls_pts = numpy.asarray(ls_pts)
    ls_bbox, ls_th, ls_ch, ls_mu, ls_el, ls_eu = read_growth_space(ls_pts)
	# Calculate error
    Ebbox = numpy.sqrt((((ls_bbox[1]-ls_bbox[0])-(targ_bbox[1]-targ_bbox[0]))/(targ_bbox[1]-targ_bbox[0]))**2 + \
			(((ls_bbox[3]-ls_bbox[2])-(targ_bbox[3]-targ_bbox[2]))/(targ_bbox[3]-targ_bbox[2]))**2) #BBox error
    Eth = numpy.sqrt(((ls_th-targ_th)/targ_th)**2) # Trunk height error
    Ech = numpy.sqrt(((ls_ch-targ_ch)/targ_ch)**2) # Crown height error
	#Eel_theta = (ls_el[0]-targ_el[0])/(2*numpy.pi)
    Eel_radii = numpy.mean(((ls_el[1:2]-targ_el[1:2])/(targ_el[1:2]))**2)
	#Eeu_theta = (ls_eu[0]-targ_eu[0])/(2*numpy.pi)
    Eeu_radii = numpy.mean(((ls_eu[1:2]-targ_eu[1:2])/(targ_eu[1:2]))**2)
    E_total = numpy.mean(numpy.array(\
		[Eth,Ech,Eel_radii,Eeu_radii]))
    return E_total

def optimise(npts,params):
# Optimise function takes in:
# npts:     number of sample points for Latin Hypercube sampling
# params:   list of parameters to be optimised (range 0:1)
# Returns variables:
# result: list of optimum values obtained from each initial condition
# error according to 'result'
    result=[];error=[];
    np = len(params)
    sample_points = pyDOE.lhs(np,samples=np*npts,criterion='center')
    for i in range(npts):
        p = sample_points[i,:]
        try:
            opt_params=optimize.minimize(estimate_error,p,method='Nelder-Mead')
            result.append(opt_params.x); error.append(opt_params.fun)
        except:
            pass
    return result,error

# Global variable
target_filename='../../obj_files/target.obj'

def main():
    print "Running main function of opt.py"

if __name__ == "__main__":
    main()

#def main():
#	result = [];
#	error = [];
#	print "Actual params were 20,5,8,75,25"
#	npts = 2000
#	sample_points = pyDOE.lhs(5, samples=npts, criterion='center')
#	for i in range(npts):
#		p = sample_points[i,:]
#		try:
#			opt_params = optimize.minimize(estimate_error,p,method='Nelder-Mead')
#			result.append(numpy.array([int(opt_params.x[0]*60), int(opt_params.x[1]*10), int(opt_params.x[2]*10), int(opt_params.x[3]*180), int(10+opt_params.x[4]*80)]))
#			error.append(opt_params.fun)
#			print numpy.array([int(opt_params.x[0]*60), int(opt_params.x[1]*10), int(opt_params.x[2]*10), int(opt_params.x[3]*180), int(10+opt_params.x[4]*80)]), opt_params.fun
#			sequence=""
#                        filename="objs/ls-%s.obj"
#                        while isfile(filename % sequence):
#                            sequence = int(sequence or 0) + 1
#                        filename = filename % sequence
#                        print filename
#                        objpts = generate_lsystem_tree_points(opt_params.x)
#                        file=open(filename,'w')
#                        for item in objpts:
#                            file.write("v %d %d %d\n" % (item[0], item[1], item[2]))
#                        file.close()
#		except:
#			pass
#
#
#	result = numpy.asarray(result)
#	error = numpy.asarray(numpy.abs(error))
#	loc = numpy.where(error == error.min())
#	print result[loc,:], error.min()
#	numpy.save('results.npy',result); numpy.save('errors.npy',error);



#if __name__ == "__main__":
#    main()
