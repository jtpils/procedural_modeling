import numpy
import pyDOE
import scipy.stats
from scipy.stats.distributions import norm
import scipy.optimize as optimize
from numpy.linalg import eig, inv
from pyquaternion import Quaternion
from os.path import isfile
import glob
#import sys
#sys.path.append('../interfacing_lysys_opt')

def generate_lsystem_tree_points(p):
	'Parse L-System string output into a list of points of cylinder base centers and their radii forming the trunk-branch skeleton representation of a tree'
	age = int(p[0]);
	no_1st_ord_branches = int(p[1]);
	no_2nd_ord_branches = int(p[2]);
	branching_angle_roll = p[3];
	branching_angle_pitch = p[4];
    #lstring = lsystem_run(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch)
    #parse l-system output string into 3d points matching Xfrog object size
    #todo

    # temporary hardcoding of tree growth: age 1-2 trunk grows taller by 2 m/yr, age 3 first branch out by 1 m/yr, age 4 second branch out by 0.5 m/yr
    # ========================
	
	trunk_growth_rate = p[5]
	branch_growth_rate = p[6]
	#trunk_growth_rate = 0.6
	#branch_growth_rate = 0.3
	pitch_angle = branching_angle_pitch/180.0 * 3.14159 #radian
	roll_angle = branching_angle_roll/180.0 * 3.14159 #radian

	current_post = numpy.array([0.0, 0.0, 0.0]) #tracking turtle's position
	current_dir = numpy.array([0.0, 0.0, 1.0])  #tracking turtle's orientation/heading (vector H)
	current_up = numpy.array([0.0, 1.0, 0.0])   #tracking turtle's up (vector U)
	current_left = numpy.array([1.0, 0.0, 0.0]) #tracking turtle's left (vector L)

	branching_stack = []

	if age == 0:
			return [current_post]

	output_point_list = [current_post]  #ground
	for t in range (1, age+1): # growing trunk
		current_post = current_post + current_dir*trunk_growth_rate   #do not use operator +=, as it overwrites all values into one current_post memory location
		output_point_list += [current_post]

	if age > 4: #1st order branches start growing
        #push state into stack
		branching_stack.append([current_post, current_dir, current_up, current_left])
        #print branching_stack

		pitch1 = Quaternion(axis=current_left,angle=pitch_angle)
		roll1 = Quaternion(axis=current_dir,angle=roll_angle)
		current_dir = pitch1.rotate(current_dir)
		current_up = pitch1.rotate(current_up)
		current_left = pitch1.rotate(current_left)
		first_branching_node = current_post

        for s in range (0, no_1st_ord_branches): #1st order branches
            current_dir = roll1.rotate(current_dir)
            current_up = roll1.rotate(current_up)
            current_left = roll1.rotate(current_left)

            for t in range(age-4): #growing 1st order branch
                current_post = first_branching_node + (t+1)*current_dir*branch_growth_rate
                output_point_list += [current_post]

            if age > 5: #2nd order branches start growing
                #push state into stack
                #print("2nd order - Current_post=%s" % (current_post))
                branching_stack.append([current_post, current_dir, current_up, current_left])
                #print branching_stack

                pitch2 = Quaternion(axis=current_left,angle=pitch_angle)
                roll2 = Quaternion(axis=current_dir,angle=roll_angle)
                current_dir = pitch2.rotate(current_dir)
                current_up = pitch2.rotate(current_up)
                current_left = pitch2.rotate(current_left)
                second_branching_node = current_post

                for u in range (0, no_2nd_ord_branches): #2nd order branch
                    current_dir = roll2.rotate(current_dir)
                    current_up = roll2.rotate(current_up)
                    current_left = roll2.rotate(current_left)

                    for v in range(age-5): #growing 2nd order branch
                        current_post = second_branching_node + (v+1)*current_dir*branch_growth_rate
                        output_point_list += [current_post]

                #pop stack to previous state
                [current_post, current_dir, current_up, current_left] = branching_stack.pop()
                #print branching_stack

        #pop stack to previous state
        [current_post, current_dir, current_up, current_left] = branching_stack.pop()
        #print branching_stack

    # ===========================================

	return output_point_list

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

def calc_growth_space_ls(pts):
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

def calc_growth_space_vx(pts):
	# Bounding box vals
	xmin,xmax = numpy.min(pts[:,0]),numpy.max(pts[:,0])
	ymin,ymax = numpy.min(pts[:,1]),numpy.max(pts[:,1])
	zmin,zmax = numpy.min(pts[:,2]),numpy.max(pts[:,2])
	bbox = numpy.asarray([xmin,xmax,ymin,xmax,zmin,zmax])
	# Estimation of trunk height
	for i in range(100):
		lh,uh = (zmax/100.)*i, (zmax/100.)*i+1
		points_at_height = numpy.sum(numpy.logical_and(pts[:,2]>lh,pts[:,2]<uh))
		if (points_at_height>250):
			trunk_height = lh;
			break;
	# Crown height
	crown_height = zmax-trunk_height;
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

#def load_target_points():
	# Function loads target_points from obj file
    #tpts = numpy.genfromtxt(target_filename,delimiter=' ',usecols=(1,2,3))
    #return tpts

def load_normalised_voxel_growth_space():
	# Function reads growth_space file, shifts x,y coordinates
	# to zero mean, sets min z coordinate to 0.
	gs = numpy.genfromtxt(gs_filename,delimiter=',')
	gs[:,0] = gs[:,0] - numpy.mean(gs[:,0])
	gs[:,1] = gs[:,1] - numpy.mean(gs[:,1])
	gs[:,2] = gs[:,2] - numpy.min(gs[:,2])
	return gs

# Global variable
target_filename='../../obj_files/target.obj'
gs_filename='../../growth-space/20171120 Tree25_VoxelCenters_10pts_25cm.csv'

tpts = load_normalised_voxel_growth_space()
targ_bbox,targ_th,targ_ch,targ_mu,targ_el,targ_eu = calc_growth_space_vx(tpts)

def estimate_error(params):
# Estimate_error function takes in:
# params: parameters used to create tree from L-System
# fname: filename of obj file (comparison data)
	# Generate L-system points for params
	ls_pts = generate_lsystem_tree_points(params)
	ls_pts = numpy.asarray(ls_pts)
	ls_bbox, ls_th, ls_ch, ls_mu, ls_el, ls_eu = calc_growth_space_ls(ls_pts)
	# Calculate error
	Ebbox = numpy.sqrt((((ls_bbox[1]-ls_bbox[0])-(targ_bbox[1]-targ_bbox[0]))/(targ_bbox[1]-targ_bbox[0]))**2 + \
		(((ls_bbox[3]-ls_bbox[2])-(targ_bbox[3]-targ_bbox[2]))/(targ_bbox[3]-targ_bbox[2]))**2) #BBox error
	Eth = numpy.sqrt(((ls_th-targ_th)/targ_th)**2) # Trunk height error
	Ech = numpy.sqrt(((ls_ch-targ_ch)/targ_ch)**2) # Crown height error
	#Eel_theta = (ls_el[0]-targ_el[0])/(2*numpy.pi)
	Eel_radii = numpy.mean(((ls_el[1:2]-targ_el[1:2])/(targ_el[1:2]))**2)
	#Eeu_theta = (ls_eu[0]-targ_eu[0])/(2*numpy.pi)
	Eeu_radii = numpy.mean(((ls_eu[1:2]-targ_eu[1:2])/(targ_eu[1:2]))**2)
	E_total = numpy.sum(numpy.array(\
		[Eth,Ech,Eel_radii,Eeu_radii]))
	#print params, E_total
	return E_total

def optimise(npts,means,stds,ranges):
# Optimise function takes in:
# npts:     number of sample points for Latin Hypercube sampling
# rflags:   array of flags to indicate whether we know ranges
# ranges:   ranges of known variables (same length as rflags)
# indices:  this should be used in a loop...
# Returns variables:
# result: list of optimum values obtained from each initial condition
# error according to 'result'
    r=[];e=[];
    nparams = len(means);
    
    # Construct sample points
    sample_points = pyDOE.lhs(nparams,samples=nparams*npts,criterion='center')
    
    # Centre around means and stds
    for i in xrange(nparams):
        sample_points[:,i] = norm(loc=means[i], scale=stds[i]).ppf(sample_points[:,i])
        
    # Test points    
    for i in range(len(sample_points)):
        p = sample_points[i,:]
		#print p
        try:
            print "Trying point:", i
            print "\tParameters:         ", p
            opt_params=optimize.minimize(estimate_error,p,method='TNC',bounds=ranges)
            print "\tOptimum parameters: ", opt_params.x
            r.append(opt_params.x); e.append(opt_params.fun)
            print "\tError:              ", opt_params.fun
        except:
            pass
    r = numpy.asarray(r); e = numpy.asarray(e);
    return r,e

def main():
# Points are in order: Ages, No 1st order, No 2nd order, Roll, pitch
        # Number of params in model
        nparams = 7;
        # Setting of default ranges - LIMITS
        default_ranges = [[0,60],[0,10],[0,10],[20,160],[0,60],[0.1,2],[0.05,3]]     
        rflags = numpy.asarray([1,1,1,0,0,0,0]) # Array to indicate number of known ranges
        ranges = default_ranges;
        # Setting of known ranges (overwrite defaults)
        ranges[0][0]=5; ranges[0][1]=30;
        ranges[1][0]=3;  ranges[1][1]=5;
        ranges[2][0]=1;  ranges[2][1]=4;
        # Setting of default means
        means = [20,4,3,90,30,0.5,0.6];
        stds = [5,1,1,20,7,0.1,0.1];
        
	print "------------------------------------------------"
	print "Running main function of opt.py..."
	print "Target values:"
	print "\tTree height:              ", targ_th+targ_ch
	print "Setting the following constraints:"
	print "\tAge:                     20-30 years"
	print "\tNo 1st order branches:   3-5"
	print "\tNo 2nd order branches:   1-4"
	print "Params with unknown range:"
	print "\tBranching angle (roll)\n\tBranching angle (pitch)";

        result,error = optimise(100,means,stds,ranges)

	result = numpy.asarray(result)
	error = numpy.asarray(numpy.abs(error))
	loc = numpy.where(error == error.min())        
                
        res = numpy.zeros(nparams)
        res[0] = result[loc,0]; res[1] = result[loc,1]
        res[2] = result[loc,2]; res[3] = result[loc,3]
        res[4] = result[loc,4]; res[5] = result[loc,5]
        res[6] = result[loc,6];

        print "Result of optimisation:"
        print "\tOptimum parameters:\t",res
        print "\tError:\t\t", error.min()
        
        output = generate_lsystem_tree_points(res)
        file = open('output.obj', 'w')
        for item in output:
            file.write("v %d %d %d\n" % (item[0], item[1], item[2]))
        file.close()

	#result = numpy.asarray(result)
	#error = numpy.asarray(numpy.abs(error))
	#loc = numpy.where(error == error.min())
	print result[loc,:], error.min()
	#numpy.save('results.npy',result); numpy.save('errors.npy',error);


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
