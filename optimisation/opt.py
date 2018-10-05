from os.path import isfile
import warnings
import random
import glob
import time
import sys

# Deal with overuse RAM
import signal,time
class TimeoutError (RuntimeError):
  pass
def handler (signum, frame):
  raise TimeoutError()
signal.signal(signal.SIGALRM,handler)
#

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
    
def write_ranges(ranges):
  
  return

# Read from xml
xml_filename='../../growth-space/example_xml/Tree1_Parameters.xml'
tpts = lgs.xml_file_gs(xml_filename,1)

def estimate_error_2(params):
  # Estimate error from rasterised tings
  # Collect points
  suppress()
  ls_mtg = interf.generate_lsystem_tree_points(params)
  allow()
  ls_pts = lgs.mtg_string_gs(ls_mtg)
  # Compute bounding boxes
  tbbox = cgs.compute_bbox(tpts)
  lsbbox = cgs.compute_bbox(ls_pts)
  # Create grids
  xmin,xmax = min(tbbox[0,0],lsbbox[0,0])-1, max(tbbox[0,1],lsbbox[0,1])
  ymin,ymax = min(tbbox[1,0],lsbbox[1,0])-1, max(tbbox[1,1],lsbbox[1,1])
  zmin,zmax = min(tbbox[2,0],lsbbox[2,0]), max(tbbox[2,1],lsbbox[2,1])
  
  xp = numpy.arange(xmin,xmax,0.5)
  yp = numpy.arange(ymin,ymax,0.5)
  zp = numpy.arange(zmin,zmax,0.5)

  t_rast  = numpy.zeros((len(xp),len(yp),len(zp)),dtype=float)
  ls_rast = numpy.zeros((len(xp),len(yp),len(zp)),dtype=float)
  
  for i in range(len(tpts)):
    x_id = numpy.abs(xp-tpts[i,0]).argmin()
    y_id = numpy.abs(yp-tpts[i,1]).argmin()
    z_id = numpy.abs(zp-tpts[i,2]).argmin()
    t_rast[x_id,y_id,z_id] = 1.0
  for j in range(len(ls_pts)):
    x_id = numpy.abs(xp-ls_pts[j,0]).argmin()
    y_id = numpy.abs(yp-ls_pts[j,1]).argmin()
    z_id = numpy.abs(zp-ls_pts[j,2]).argmin()
    ls_rast[x_id,y_id,z_id] = 1.0
  
  E = numpy.sum(numpy.sum(numpy.sum(numpy.abs(t_rast-ls_rast))))/numpy.sum(numpy.sum(numpy.sum(t_rast)))
  
  comments = '# '
  
  if (numpy.abs((lsbbox[2,1]-tbbox[2,1]))/tbbox[2,1])>0.2:
    E = E+1.0
    comments = comments + "tree differs in height |"
    
  if (lsbbox[2,1]<0.0):
    comments = comments + " branches intersect ground |"
    E = E+1.0
    
  return E,comments

def create_sample_points(npts,means,stds,ranges):
	# Number of parameters
	nparams = len(means);
	# Construct initial sample points
	temp_points = pyDOE.lhs(nparams,samples=nparams*npts)
	# Centre around means and stds
	for i in range(len(temp_points)):
	  for j in xrange(nparams):
	    temp_points[i,j] = ranges[j,0] + temp_points[i,j]*(ranges[j,1]-ranges[j,0])
	sample_points = []
	for i in range(len(temp_points)):
	  add_flag=0
	  for j in xrange(nparams):
	    if temp_points[i,j]>ranges[j,1] or temp_points[i,j]<ranges[j,0]:
	      print j, "issue"
	      add_flag=add_flag+1
	  if add_flag==0:
	    sample_points.append(temp_points[i,:])
	
	sample_points = numpy.asarray(sample_points)
	return sample_points

def map_error(points,fname):
	p=[];e=[];
	for i in range(len(points)):
	  stat_str = "Testing point %i/%i" % (i,len(points))
	  sys.stdout.write('%s\r' % stat_str)
	  sys.stdout.flush()
	  try:
	    signal.alarm(20)
    	    #print "\tPoint: ", points[i,:]
	    err,comments = estimate_error_2(points[i,:])
	    e.append(err)
	    p.append(points[i,:])
	    #print "\tError: ", err 
	    signal.alarm(0)
	    pstr = '%i\t%0.2f\t%0.2f\t%2.2f\t%i\t%i\t%2.2f\t%2.2f\t%0.2f\t%i\t%0.2f\t%f\t%s\n' % \
	      (points[i,0], points[i,1], points[i,2], points[i,3], points[i,4], points[i,5], \
		points[i,6], points[i,7], points[i,8], points[i,9], points[i,10], err, comments)
	    with open(fname,'a+') as f:
	      f.write(pstr)
	  except TimeoutError as ex:
	    #print "Point timed out"
	    continue
	return p,e    
    
def select_from_population(population, population_error, num_best, num_lucky):
    sort_indices = numpy.argsort(population_error);
    sorted_population = population[sort_indices,:]
    sorted_error = population_error[sort_indices]
    next_sample_points = [];
    chosen_errors = [];
    for i in range(num_best):
      next_sample_points.append(sorted_population[i,:]);
      chosen_errors.append(sorted_error[i])
    for i in range(num_lucky):
      ind = random.randint(0,len(sorted_population)-1)
      next_sample_points.append(sorted_population[ind,:]);
      chosen_errors.append(sorted_error[ind])
    return next_sample_points, chosen_errors
  
def create_children(parent1,parent2,inheritance):
    child1 = parent1; child2 = parent2;
    child1[inheritance] = parent2[inheritance];
    child2[inheritance] = parent1[inheritance];
    return child1, child2
  
def population_breeding(population,numchildren):
    next_population = []
    for i in range(int(numpy.floor(len(population)/2))):
      for j in range(numchildren):
	parent1 = population[i]; parent2 = population[len(population)-1-i]
	inheritance = numpy.random.choice([0,1], size=(len(parent1)), p=[0.5,0.5])
	child1, child2 = create_children(parent1,parent2,inheritance)
	next_population.append(child1)
	next_population.append(child2)
    return next_population
  
def mutate_child(child,amp):
    mutation_flag = numpy.random.choice([-1,1], size=(len(child)), p=[0.5,0.5])
    for i in range(len(child)):
	if i==0 or i==4 or i==5:
	  if child[i]<=1:
	    child[i]=child[i]+numpy.abs(mutation_flag[i])
	  else:
	    child[i]=child[i]+mutation_flag[i]
	else:
	    child[i]=child[i]*(1+amp*mutation_flag[i])
    return child
  
def mutate_population(population,amp,chance):
    mutated_population=[]
    for i in range(len(population)):
      if random.random()<chance:
	mutated_population.append(mutate_child(population[i],amp))
    return mutated_population

def main():
# Points are in order: Ages, No 1st order, No 2nd order, Roll, pitch
    # Number of params in model
    params = ['Age', 'Trunk pitch angle', 'Trunk roll angle', 'Trunk height',\
      'No. 1st order branches', 'No. 2nd order branches', 'Branch pitch angle',\
	'Branch roll angle', 'Diameter growth rate', 'Annual no. new nodes',\
	  'Average internode length']

    default_ranges = [\
        [5,30], # age = params[0
        [-10,10],[-10,10],[1,7],[1,4],[1,4],\
        [10,70],[30,190],[0.01,1],[5,60],[0.01,0.1]]
    
    nparams = len(default_ranges);
    ranges = default_ranges;
    
    means = [20,0,0,6,2,2,30,30,0.05,30,0.03]
    stds = [7,0.5,0.5,2,1,1,1,15,15,0.025,10,0.01];

    print "------------------------------------------------"
    print "Running main function of opt.py..."
    print "\tTarget values set from file   : %s" % xml_filename
    print "Setting the following ranges:"
    for i in range(len(params)):
      print "\t{:<30}:\t{:>6.2f} - {:<6.2f}".format(params[i], ranges[i][0], ranges[i][1])
    print "------------------------------------------------"

    ranges = numpy.asarray(ranges)

    run=2; save_obj=1;
    
    npoints = int(raw_input("Select number of sample points per variable for initial population: "));
    ef_name = 'func.dat'
    
    if run==2:
	# Create initial population
        initial_sample_points = create_sample_points(npoints,means,stds,ranges)        
        print "Testing initial %i sample points" % (npoints*nparams)	
	print "Outputting results to file: %s" % ef_name
	res = []; err=[];
	t0 = time.time()
	try:
	  rpts,error=map_error(initial_sample_points,ef_name);
	except:
	  pass
	# Load all historic data from file
	t1 = time.time()
	print "Initial %i points took %d seconds" % (npoints*nparams,(t1-t0))
	print "------------------------------------------------"
	results = numpy.loadtxt(ef_name,delimiter='\t',usecols=(0,1,2,3,4,5,6,7,8,9,10,11))
	error = results[:,-1]
	results = results[:,0:11]
	total_pop_size = len(error);
	# Params for selecting generations
	num_generations = 50#int(raw_input("Select number of generations: "))
	numbest = 75; numrandom=25;
	# Params for mutations
	amp=0.1; mut_chance=0.75;
	for i in range(num_generations):
	  t2 = time.time()
	  print "Performing generation number: %i" % (i+1)
	  results = numpy.loadtxt(ef_name,delimiter='\t',usecols=(0,1,2,3,4,5,6,7,8,9,10,11))
	  error = results[:,-1]; results = results[:,0:11]; total_pop_size = len(error);
	  print "\tSelecting %i best performing and %i random individuals from population of %i" % (numbest, numrandom, total_pop_size)
	  parent_sample_points, parent_errors = select_from_population(results,error,numbest,numrandom)
  	  print "\tMinimum error in parent population is %f" % min(parent_errors)
	  print "\tMaximum error in parent population is %f" % max(parent_errors)
	  child_sample_points = population_breeding(parent_sample_points,2);
	  mutated_sample_points = mutate_population(child_sample_points,amp,mut_chance)
	  print "\tTesting %i sample points from new population" % len(mutated_sample_points)
	  #try:
	  rpts,error = map_error(numpy.asarray(mutated_sample_points),ef_name)
	  #except:
	  #  continue
	  t3 = time.time();
	  print "\tTesting generation %i points took %d seconds" % (i+1,(t3-t2))

	print "------------------------------------------------"
	print "Optimisation complete in %d seconds" % ((t3-t0)/60);
	print "------------------------------------------------"
      
if __name__ == "__main__":
    main()
    
    