import numpy
import pyDOE
import scipy.stats
from scipy.stats.distributions import norm
import scipy.optimize as optimize
from pyquaternion import Quaternion
from os.path import isfile
import glob
import multiprocessing as mp
import sys
import warnings
warnings.filterwarnings("ignore")
import calc_growth_space as cgs
import load_growth_space as lgs

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

# Global variable
#target_filename='../../obj_files/target.obj'
#gs_filename='../../growth-space/20171120 Tree25_VoxelCenters_10pts_25cm.csv'
gs_dir = '../../growth-space/voxel_size_tests/'
gs_filename='Tree1_50cm_5pts.csv'
resolution=50

# Uncomment to use obj as target
# Actual params were 20,5,8,75,25,2,1
# tpts = numpy.genfromtxt(target_filename,delimiter=' ',usecols=(1,2,3))
# targ_bbox,targ_th,targ_ch,targ_mu,targ_el,targ_eu = calc_growth_space_ls(tpts)

tpts = lgs.normalised_voxel_gs(gs_dir+gs_filename,resolution)
targ_bbox,targ_th,targ_ch,targ_mu,targ_el,targ_eu = cgs.cgs_vx(tpts)

def estimate_error_1(params):
# Estimate_error function takes in:
# params: parameters used to create tree from L-System
# fname: filename of obj file (comparison data)
    # Generate L-system points for params
    ls_pts = generate_lsystem_tree_points(params)
    ls_pts = numpy.asarray(ls_pts)
    ls_bbox, ls_th, ls_ch, ls_mu, ls_el, ls_eu = cgs.cgs_ls(ls_pts)
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
# Estimate_error function takes in:
# params: parameters used to create tree from L-System
# fname: filename of obj file (comparison data)
    # Generate L-system points for params
    ls_pts = generate_lsystem_tree_points(params)
    ls_pts = numpy.asarray(ls_pts)
    ls_bbox, ls_th, ls_ch, ls_mu, ls_el, ls_eu = cgs.cgs_ls(ls_pts)
    # Calculate error
    Ebbox = numpy.sqrt((((ls_bbox[1]-ls_bbox[0])-(targ_bbox[1]-targ_bbox[0]))/(targ_bbox[1]-targ_bbox[0]))**2 + \
        (((ls_bbox[3]-ls_bbox[2])-(targ_bbox[3]-targ_bbox[2]))/(targ_bbox[3]-targ_bbox[2]))**2) #BBox error
    Eth = numpy.sqrt(((ls_th-targ_th)/targ_th)**2) # Trunk height error
    Ech = numpy.sqrt(((ls_ch-targ_ch)/targ_ch)**2) # Crown height error
    #Eel_theta = (ls_el[0]-targ_el[0])/(2*numpy.pi)
    Eel_radii = numpy.mean(((ls_el[1:2]-targ_el[1:2])/(targ_el[1:2]))**2)
    #Eeu_theta = (ls_eu[0]-targ_eu[0])/(2*numpy.pi)
    Eeu_radii = numpy.mean(((ls_eu[1:2]-targ_eu[1:2])/(targ_eu[1:2]))**2)

    E_mu = numpy.mean((ls_mu-targ_mu)**2)

    E_total = numpy.sum(numpy.array(\
        [Eth,Ech,Eel_radii,Eeu_radii,E_mu]))
    #print params, E_total
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
#		opt_params=optimize.minimize(estimate_error_1,points,\
#				method='SLSQP',bounds=ranges)

		opt_params=optimize.minimize(estimate_error_1,points,\
				method='TNC',bounds=ranges)
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
# Points are in order: Ages, No 1st order, No 2nd order, Roll, pitch
	# Number of params in model
	nparams = 7;
	# Setting of default ranges - LIMITS
	default_ranges = [[0,60],[0,10],[0,10],[20,160],[0,60],[0.1,3],[0.05,3]]
	rflags = numpy.asarray([1,1,1,0,0,0,0]) # Array to indicate number of known ranges
	ranges = default_ranges;
	# Setting of known ranges (overwrite defaults)
	ranges[0][0]=5; ranges[0][1]=30;
	ranges[1][0]=2;  ranges[1][1]=6;
	ranges[2][0]=1;  ranges[2][1]=8;
	# Setting of default means
	means = [20,3,3,75,25,2.,1.];
	stds = [3,1,1,25,10,0.5,0.25];

	print "------------------------------------------------"
	print "Running main function of opt.py..."
	print "Target values:"
	print "\tTrunk height:              ", targ_th
	print "\tCrown height:              ", targ_ch
	print "Setting the following constraints:"
	print "\tAge:                 :   5-30"
	print "\tNo 1st order branches:   2-4"
	print "\tNo 2nd order branches:   1-4"
	print "Params with unknown range:"
	print "\tBranching angle (roll)\n\tBranching angle (pitch)";

	run=1; save_npy=0; save_obj=1;
	npoints = int(sys.argv[1]);
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
			out_file = 'opt_npts' + str(npoints) + '.obj'
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
