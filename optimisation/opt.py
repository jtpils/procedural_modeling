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
#sys.path.append('../interfacing_lysys_opt')

def generate_lsystem_tree_points(p):
	'Parse L-System string output into a list of points of cylinder base centers and their radii forming the trunk-branch skeleton representation of a tree'
	age = int(p[0]);
	no_1st_ord_branches = int(p[1]);
	no_2nd_ord_branches = int(p[2]);
	branching_angle_roll = p[3];
	branching_angle_pitch = p[4];
	trunk_growth_rate = p[5]
	branch_growth_rate = p[6]
    # temporary hardcoding of tree growth: age 1-2 trunk grows taller by 2 m/yr, age 3 first branch out by 1 m/yr, age 4 second branch out by 0.5 m/yr
    # ========================
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

def load_mtg_file(fname):
	mtg_file = fname
	print "Reading file: %s" % mtg_file
	with open(mtg_file) as f:
		lines = f.readlines()
	i=0
	while(1):
		if lines[i].startswith("ENTITY-CODE"):
			break
		else:
			i=i+1
	# Columns are 1=y, 2=radius, 3=z, 4=z
	pts = np.genfromtxt(mtg_file,skip_header=i+1,usecols=(4,1,3))
	return pts

def load_normalised_voxel_growth_space():
	# Function reads growth_space file, shifts x,y coordinates
    # to zero mean, sets min z coordinate to 0.
    gs = numpy.genfromtxt(gs_dir+gs_filename,delimiter=',')
    gs = gs*resolution/100.
    gs[:,0] = gs[:,0] - numpy.mean(gs[:,0])
    gs[:,1] = gs[:,1] - numpy.mean(gs[:,1])
    gs[:,2] = gs[:,2] - numpy.min(gs[:,2])
    return gs

# Global variable
#target_filename='../../obj_files/target.obj'
#gs_filename='../../growth-space/20171120 Tree25_VoxelCenters_10pts_25cm.csv'
gs_dir = '../../growth-space/voxel_size_tests/'
gs_filename=sys.argv[1]
resolution=int(sys.argv[2])


# Uncomment to use obj as target
# Actual params were 20,5,8,75,25,2,1
# tpts = numpy.genfromtxt(target_filename,delimiter=' ',usecols=(1,2,3))
# targ_bbox,targ_th,targ_ch,targ_mu,targ_el,targ_eu = calc_growth_space_ls(tpts)

tpts = load_normalised_voxel_growth_space()
targ_bbox,targ_th,targ_ch,targ_mu,targ_el,targ_eu = calc_growth_space_vx(tpts)

def estimate_error_1(params):
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

def estimate_error_2(params):
# To do

    return E_total

def create_sample_points(npts,means,stds):
	# Number of parameters
	nparams = len(means);
	# Construct initial sample points
	sample_points = pyDOE.lhs(nparams,samples=nparams*npts,criterion='center')
	# Centre around means and stds
	for i in xrange(nparams):
		sample_points[:,i] = norm(loc=means[i], scale=stds[i]).ppf(sample_points[:,i])
	return sample_points

def optimise(points,ranges):
	r=[]; e=[]; i=[];
	try:
		opt_params=optimize.minimize(estimate_error_1,points,\
				method='SLSQP',bounds=ranges)
		#print "\tOpt params: ", opt_params.x
		r = opt_params.x; e = opt_params.fun
		#print "\tError:      ", opt_params.fun
		i = points
	except:
		r = points
		e = 100.0;
		i = points
		pass
	i = numpy.asarray(i)
	r = numpy.asarray(r); e = numpy.asarray(e);
	return i,r,e

def main():
        print sys.argv[1]
# Points are in order: Ages, No 1st order, No 2nd order, Roll, pitch
	# Number of params in model
	nparams = 7;
	# Setting of default ranges - LIMITS
	default_ranges = [[0,60],[0,10],[0,10],[20,160],[0,60],[0.1,3],[0.05,3]]
	rflags = numpy.asarray([1,1,1,0,0,0,0]) # Array to indicate number of known ranges
	ranges = default_ranges;
	# Setting of known ranges (overwrite defaults)
	ranges[0][0]=5; ranges[0][1]=30;
	ranges[1][0]=3;  ranges[1][1]=5;
	ranges[2][0]=1;  ranges[2][1]=10;
	# Setting of default means
	means = [20,5,8,75,25,2.,1.];
	stds = [3,2,2,25,10,0.5,0.25];

	print "------------------------------------------------"
	print "Running main function of opt.py..."
	print "Target values:"
	print "\tTree height:              ", targ_th+targ_ch
	print "Setting the following constraints:"
	print "\tAge:                     5-30 years"
	print "\tNo 1st order branches:   3-5"
	print "\tNo 2nd order branches:   1-10"
	print "Params with unknown range:"
	print "\tBranching angle (roll)\n\tBranching angle (pitch)";

	run=1; save_npy=0; save_obj=1;
	npoints = 100;
	num_threads=mp.cpu_count();


	if run==1:
		# Run optimisation routine
		sample_points = create_sample_points(npoints,means,stds)
		print "Testing %i sample points with %i cpus" % (len(sample_points), num_threads)
		# Create pool of threads
		pool = mp.Pool(processes=num_threads);
		# Run on nprocs
		result = [ pool.apply_async(optimise,args=(i,ranges,)) for i in sample_points]
		# Collect data when ready
		output = numpy.asarray([p.get() for p in result])
		# Reorganise results
		initial_points = numpy.asarray(output[:,0])
		results = numpy.asarray(output[:,1])
		error = numpy.asarray(output[:,2])
		loc = numpy.where(error==error.min())
		loc = int(loc[0])
		#Convert result to array (don't remember why it was necessary
		# to do it like this
		res = numpy.zeros(nparams)
		res[0] = results[loc][0]; res[1] = results[loc][1]
		res[2] = results[loc][2]; res[3] = results[loc][3]
		res[4] = results[loc][4]; res[5] = results[loc][5]
		res[6] = results[loc][6];
		# Print results to screen
		print "Result of optimisation:"
		print "\tOpt params:",res
		print "\tError:\t", error.min()
		# Save as obj?
		if save_obj==1:
			output = generate_lsystem_tree_points(res)
			out_file = 'opt_' + gs_filename + '.obj'
			file = open(out_file, 'w')
			for item in output:
				file.write("v %d %d %d\n" % (item[0], item[1], item[2]))
			file.close()
		# Save as npy
		if save_npy==1:
			numpy.save('results.npy',result); numpy.save('errors.npy',error);
			numpy.save('initial_points.npy',initials)
	print "------------------------------------------------"


if __name__ == "__main__":
    main()
