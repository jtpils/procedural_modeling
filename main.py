from os.path import isfile
import warnings
import signal
import time
import numpy
import pyDOE
import sys

#import optimisation.calc_growth_space as cgs
import optimisation.load_growth_space as lgs
import optimisation.read_xml as rxml
import optimisation.opt as opt

warnings.filterwarnings("ignore")

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

def write_optimums(fname,opt_params,params,err):
  f = open(fname,"w")
  hstr = '# Optimum parameter configuration: E=%f\n' % err
  f.write(hstr)
  for i in range(len(params)):
    ostr = '{:<6.2f}\t# {}\n'.format(opt_params[i],params[i])
    f.write(ostr)
  return

def main():
    params = ['Age', 'Trunk pitch angle', 'Trunk roll angle', 'Trunk height',\
      'No. 1st order branches', 'No. 2nd order branches', 'Branch pitch angle',\
	'Branch roll angle', 'Diameter growth rate', 'Annual no. new nodes',\
	  'Average internode length']

    default_ranges = [\
        [5,20], # age = params[0
        [-5,5],[-5,5],[2,6],[2,4],[2,4],\
        [20,40],[60,200],[0.1,1.0],[10,20],[0.1,2]]

    nparams = len(default_ranges);
    ranges = default_ranges;

    means = [20,0,0,6,2,2,30,30,0.05,30,0.03]
    stds = [7,0.5,0.5,2,1,1,1,15,15,0.025,10,0.01];

    # Read from xml
    xml_filename='../growth-space/example_xml/Tree1_Parameters.xml'
    tpts = lgs.xml_file_gs(xml_filename,1)

    print "------------------------------------------------"
    print "Running main function of opt.py..."
    print "\tTarget values set from file   : %s" % xml_filename
    cname, species, location, th = rxml.rxml_treeparams(xml_filename)
    print "\tTree common name : %s" % cname
    print "\tTree species     : %s" % species
    print "\tTree height      : %d m" % th
    print "Setting the following ranges:"
    for i in range(len(params)):
      print "\t{:<30}:\t{:>6.2f} - {:<6.2f}".format(params[i], ranges[i][0], ranges[i][1])
    print "------------------------------------------------"

    ranges = numpy.asarray(ranges)

    run=1; save_obj=1;

    npoints = int(raw_input("Select number of sample points per variable for initial population: "));
    ef_name = 'func.dat'

    if run==1:
      	t0 = time.time()
	if npoints>0:
	  # Create initial population
	  initial_sample_points = opt.create_sample_points(npoints,means,stds,ranges)
	  print "Testing initial %i sample points" % len(initial_sample_points)
	  print "Outputting results to file: %s" % ef_name
	  res = []; err=[];
	  try:
	    rpts,error=opt.map_error(initial_sample_points,tpts,ef_name);
	  except:
	    pass
	  # Load all historic data from file
	  t1 = time.time()
	  print "Initial %i points took %d seconds" % (npoints*nparams,(t1-t0))
	  print "------------------------------------------------"
	else:
	  print "Skipping initial sample points and loading from results file: %s" % ef_name
	results = numpy.loadtxt(ef_name,delimiter='\t',usecols=(0,1,2,3,4,5,6,7,8,9,10,11))
	error = results[:,-1]
	results = results[:,0:11]
	total_pop_size = len(error);
	# Params for selecting generations
	num_generations = 50#int(raw_input("Select number of generations: "))
	numbest = 90; numrandom=10;
	# Params for mutations
	amp=0.5; mut_chance=1.0;
	weights = [5,4,3,2,1];
	for i in range(num_generations):
	  t2 = time.time()
	  print "Performing generation number: %i/%i" % ((i+1),num_generations)
	  results = numpy.loadtxt(ef_name,delimiter='\t',usecols=(0,1,2,3,4,5,6,7,8,9,10,11))
	  error = results[:,-1]; results = results[:,0:11]; total_pop_size = len(error);
	  print "\tSelecting %i best performing and %i random individuals from population of %i" % (numbest, numrandom, total_pop_size)
	  parent_sample_points, parent_errors = opt.select_from_population(results,error,numbest,numrandom)
	  write_optimums("opt_params.txt",parent_sample_points[0],params,parent_errors[0])
  	  print "\tMinimum error in parent population is %f" % min(parent_errors)
	  print "\tMaximum error in parent population is %f" % max(parent_errors)
	  clean_parent_sample_points = opt.check_population(parent_sample_points,ranges)
	  extended_points = opt.extend_population(clean_parent_sample_points,weights)
	  child_sample_points = opt.population_breeding(extended_points,2)
	  clean_child_sample_points = opt.check_population(child_sample_points,ranges)
	  mutated_sample_points = opt.mutate_population(clean_child_sample_points,amp,mut_chance)
	  clean_mutated_sample_points = opt.check_population(mutated_sample_points,ranges)
	  print "\tTesting %i sample points from new population" % len(clean_mutated_sample_points)
	  try:
	    rpts,error = opt.map_error(numpy.asarray(clean_mutated_sample_points),tpts,ef_name)
	  except:
	    continue
	  t3 = time.time();
	  print "\tTesting generation %i points took %d seconds" % (i+1,(t3-t2))


	results = numpy.loadtxt(ef_name,delimiter='\t',usecols=(0,1,2,3,4,5,6,7,8,9,10,11))
	error = results[:,-1]; results = results[:,0:11]; total_pop_size = len(error);
	parent_sample_points, parent_errors = opt.select_from_population(results,error,5,0)
	write_optimums("opt_params.txt",parent_sample_points[0],params,parent_errors[0])
	print "------------------------------------------------"
	print "Optimisation complete in %d minutes" % ((t3-t0)/60);
	print "Optimum parameter configuration from %i tested parameter combinations is:" % len(error)
	for i in range(len(params)):
	  print "\t{:<30}:\t{:<6.2f}".format(params[i], parent_sample_points[0][i])
	print "------------------------------------------------"


if __name__ == "__main__":
    main()