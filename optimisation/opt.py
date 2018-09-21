from os.path import isfile
import warnings
import glob
import sys
warnings.filterwarnings("ignore")
sys.path.append('../interfacing_lsys_opt/')


import calc_growth_space as cgs
import load_growth_space as lgs
import lsystem_temp as lsys
import interface as interf

import multiprocessing as mp
from scipy.stats.distributions import norm
import scipy.optimize as optimize
import scipy.stats
import numpy
import pyDOE

class NullWriter(object):
    def write(self,arg):
        pass
    def flush(self):
        pass
nullwrite = NullWriter()
oldstdout = sys.stdout

def suppress():
    sys.stdout = nullwrite

def allow():
    sys.stdout = oldstdout


# Global variable
gs_dir = '../../growth-space/voxel_size_tests/'
gs_filename='Tree1_50cm_5pts.csv'
resolution=50
xml_filename='../../growth-space/example_xml/Tree1_Parameters.xml'

# Read from gs csv
#tpts = lgs.normalised_voxel_gs(gs_dir+gs_filename,resolution)
#print numpy.shape(tpts)

# Read from xml
tpts = lgs.xml_file_gs(xml_filename)
#print numpy.shape(tpts)

targ_bbox,targ_th,targ_ch,targ_mu,targ_el,targ_eu = cgs.cgs_vx(tpts)


def estimate_error_1(params):
# Estimate_error function takes in:
# params: parameters used to create tree from L-System
# fname: filename of obj file (comparison data)
    # Generate L-system points for params
    suppress()
    ls_pts = interf.generate_lsystem_tree_points(params)
    allow()
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

def create_sample_points(npts,means,stds):
	# Number of parameters
	nparams = len(means);
	# Construct initial sample points
	sample_points = pyDOE.lhs(nparams,samples=nparams*npts,criterion='center')
	# Centre around means and stds
	for i in xrange(nparams):
		sample_points[:,i] = norm(loc=means[i], scale=stds[i]).ppf(sample_points[:,i])
	return sample_points

def optimise_p(points,ranges,method):
	r=[]; e=[]; i=[];
	try:
		suppress()
		if method=='SLSQP':
			opt_params=optimize.minimize(estimate_error_1,points,\
				method='SLSQP',bounds=ranges)
		elif method=='TNC':
			opt_params=optimize.minimize(estimate_error_1,points,\
				method='TNC',bounds=ranges)
		else:
			allow()
			print "Optimisation method not recognised."
		allow()
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

def optimise_s(points,ranges,method):
	r=[]; e=[]; i=[];
	for j in range(len(points)):
		try:
			suppress()
			if method=='SLSQP':
				opt_params=optimize.minimize(estimate_error_1,points[j,:],\
				    method='SLSQP',bounds=ranges)
			elif method=='TNC':
				opt_params=optimize.minimize(estimate_error_1,points[j,:],\
				    method='TNC',bounds=ranges)
			else:
				print "Optimisation method not recognised."
		        allow()
			print "\tPoint ", j, ": Opt params: ", opt_params.x
			r.append(opt_params.x); e.append(opt_params.fun)
			print "\t\tError:      ", opt_params.fun
			i.append(points[j,:])
		except:
		        print "Not working"
			r.append(points[j,:])
			e.append(100.0)
			i.append(points[j,:])
			pass
	i = numpy.asarray(i);
	r = numpy.asarray(r);
	e = numpy.asarray(e);		      
	return i,r,e

def main():
# Points are in order: Ages, No 1st order, No 2nd order, Roll, pitch
    # Number of params in model
    
    # age = params[0]
    # trunk_pitch_angle=params[1]
    # trunk_roll_angle=params[2],
    # trunk_height=params[3],
    # no_first_ord_branches=int(params[4]),
    # no_second_ord_branches=int(params[5]),
    # branching_pitch_angle=params[6],
    # branching_roll_angle=params[7],
    # diameter_growth_rate=params[8],
    # annual_no_new_nodes=params[9],
    # avg_internode_length=params[10]
	# Setting of default ranges - LIMITS
    default_ranges = [\
        [0,60], # age = params[0
        [-10,10],[-10,10],[1,20],[1,4],[1,4],\
        [10,70],[10,70],[0.01,1],[10,50],[0.01,0.1]]
    nparams = len(default_ranges);
    ranges = default_ranges;
    means = [20,0,0,6,2,2,30,30,0.05,30,0.03]
	# Setting of known ranges (overwrite defaults)
	# ranges[0][0]=5; ranges[0][1]=30;
	# ranges[1][0]=2;  ranges[1][1]=6;
	# ranges[2][0]=1;  ranges[2][1]=8;
	# Setting of default means
    stds = [7,0.5,0.5,2,1,1,1,15,15,0.025,10,0.01];

    print "------------------------------------------------"
    print "Running main function of opt.py..."
    print "Target values:"
    print "\tTrunk height:              ", targ_th
    print "\tCrown height:              ", targ_ch
    print "Setting the following constraints:"
    print "\tAge:                       ", ranges[0][0],"-",ranges[0][1]
    print "\tTrunk pitch angle:         ", ranges[1][0],"-",ranges[1][1]
    print "\tTrunk roll angle:          ", ranges[2][0],"-",ranges[2][1]
    print "\tTrunk height:              ", ranges[3][0],"-",ranges[3][1]
    print "\tNo. 1st order branches:    ", ranges[4][0],"-",ranges[4][1]
    print "\tNo. 2nd order branches:    ", ranges[5][0],"-",ranges[5][1]
    print "\tBranch pitch angle:        ", ranges[6][0],"-",ranges[6][1]
    print "\tBranch roll angle:         ", ranges[7][0],"-",ranges[7][1]
    print "\tDiameter growth rate:      ", ranges[8][0],"-",ranges[8][1]
    print "\tAnnual no. new nodes:      ", ranges[9][0],"-",ranges[9][1]
    print "\tAverage internode length:  ", ranges[10][0],"-",ranges[10][1]
    print "------------------------------------------------"
    
    run=2; save_npy=0; save_obj=1;
    method='Buckshot';#'SLSQP'
    npoints = int(raw_input("Select number of sample points per variable: "));
    avail_threads=mp.cpu_count();
    thread_string = "Select number of processors (%i available): " % avail_threads
    num_threads = int(raw_input(thread_string))

    if run==1:
		# Run optimisation routine
        sample_points = create_sample_points(npoints,means,stds)
        print "Testing %i sample points with %i cpus" % (len(sample_points), num_threads)
		# Create pool of thread
        if num_threads>1:
	    print "Running in parallel"
	    print "------------------------------------------------"
            pool = mp.Pool(processes=num_threads,maxtasksperchild=1000);
		    # Run on nprocs
            result = [ pool.apply_async(optimise_p,args=(i,ranges,method)) for i in sample_points]
            output = numpy.asarray([p.get() for p in result])
            initial_points = numpy.asarray(output[:,0])
            results = numpy.asarray(output[:,1])
            error = numpy.asarray(output[:,2])
        else:
	    print "Running in serial"
	    print "------------------------------------------------"
            initial_points, results, error = optimise_s(sample_points,ranges,method)
        allow()
        # Collect data when ready

		# Reorganise results

        loc = numpy.where(error==error.min())
        loc = int(loc[0])
		#Convert result to array (don't remember why it was necessary
		# to do it like this
        res = numpy.zeros(nparams)
        for i in range(nparams):
		res[i] = results[loc][i];
		# Print results to screen
        print "Result of optimisation:"
        print "\tOpt params:",res
        print "\tError:\t", error.min()
		# Save as obj?
        if save_obj==1:
			output = interf.generate_lsystem_tree_points(res)
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
    
    if run==2:
        sample_points = create_sample_points(npoints,means,stds)
        from mystic.solvers import BuckshotSolver
        from mystic.solvers import PowellDirectionalSolver
        from mystic.termination import NormalizedChangeOverGeneration as NCOG
	try:
	  from pathos.pools import ProcessPool as Pool
	except ImportError:
	  from mystic.pools import SerialPool as Pool
	solver = BuckshotSolver(len(ranges),npoints)
	solver.SetNestedSolver(PowellDirectionalSolver)
	solver.SetMapper(Pool().map)
	solver.SetInitialPoints(sample_points)
)
	solver.Solve(estimate_error_1, NCOG(1e-4), disp=1)
	suppress()
	solution = solver.Solution()
	allow()



if __name__ == "__main__":
    main()
