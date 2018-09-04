from os.path import isfile
import warnings
import glob
import sys
warnings.filterwarnings("ignore")

import calc_growth_space as cgs
import load_growth_space as lgs
import lsystem_temp as lsys

import multiprocessing as mp
from scipy.stats.distributions import norm
import scipy.optimize as optimize
import scipy.stats
import numpy
import pyDOE

# Global variable
gs_dir = '../../growth-space/voxel_size_tests/'
gs_filename='Tree1_50cm_5pts.csv'
resolution=50

tpts = lgs.normalised_voxel_gs(gs_dir+gs_filename,resolution)
targ_bbox,targ_th,targ_ch,targ_mu,targ_el,targ_eu = cgs.cgs_vx(tpts)

def estimate_error_1(params):
# Estimate_error function takes in:
# params: parameters used to create tree from L-System
# fname: filename of obj file (comparison data)
    # Generate L-system points for params
    ls_pts = lsys.generate_lsystem_tree_points(params)
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
    ls_pts = lsys.generate_lsystem_tree_points(params)
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

def optimise(points,ranges,method):
	r=[]; e=[]; i=[];
	try:
		if method=='SLSQP':
			opt_params=optimize.minimize(estimate_error_1,points,\
				method='SLSQP',bounds=ranges)
		elif method=='TNC':
			opt_params=optimize.minimize(estimate_error_1,points,\
				method='TNC',bounds=ranges)
		else:
			print "Optimisation method not recognised."
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

	run=1; save_npy=0; save_obj=1; method='SLSQP'
	npoints = int(sys.argv[1]);
	num_threads=mp.cpu_count();


	if run==1:
		# Run optimisation routine
		sample_points = create_sample_points(npoints,means,stds)
		print "Testing %i sample points with %i cpus" % (len(sample_points), num_threads)
		# Create pool of threads
		pool = mp.Pool(processes=num_threads);
		# Run on nprocs
		result = [ pool.apply_async(optimise,args=(i,ranges,method)) for i in sample_points]
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
			output = lsys.generate_lsystem_tree_points(res)
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
