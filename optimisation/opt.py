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
#from interfacing_lsys_opt import interface as interf
import interfacing_lsys_opt.interface as interf

import calc_growth_space as cgs
import load_growth_space as lgs
import read_xml as rxml

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


def estimate_error_2(sp_id,params,tpts):
  # Estimate error from rasterised tings
  # Collect points
  suppress()
  ls_mtg = interf.generate_lsystem_tree_points(sp_id,params)
  allow()
  ls_pts = lgs.mtg_string_gs(ls_mtg)/10.
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

  locs = numpy.logical_and(t_rast,ls_rast)
# EF1
  E1 = 1-numpy.sum(numpy.sum(numpy.sum(numpy.logical_and(t_rast,ls_rast))))/numpy.sum(numpy.sum(numpy.sum(t_rast)))
# EF2
  E2 = numpy.sum(numpy.sum(numpy.sum(numpy.abs(t_rast-ls_rast))))/numpy.sum(numpy.sum(numpy.sum(t_rast)))
  comments = '# '

  E3 = 0.0

  if (numpy.abs((lsbbox[0,1]-lsbbox[0,0])-(tbbox[0,1]-tbbox[0,0]))/(tbbox[0,1]-tbbox[0,0])>0.5) \
    or (numpy.abs((lsbbox[1,1]-lsbbox[1,0])-(tbbox[1,1]-tbbox[1,0]))/(tbbox[1,1]-tbbox[1,0])>0.5):
    E3 = E3 + 1.0
    comments = comments + "LS/GS width difference | "

  if numpy.abs(lsbbox[2,1]-tbbox[2,1])/tbbox[2,1]>0.25:
    E3 = E3 + 1.0
    app_str = "LS/GS height difference | "
    comments = comments + app_str

  if (lsbbox[2,0]<0):
    E3 = E3 + 1.0
    comments = comments + "LS tree intersects ground |"

  app_str = "LS tree height=%f | GS tree height=%f" % (lsbbox[2,1], tbbox[2,1])
  comments = comments + app_str



  E = E1 + 0.1*E2 + E3

  return E, comments

def create_sample_points(npts,ranges):
	# Number of parameters
	nparams = len(ranges);
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

def map_error(species_id,points,tpoints,fname):
	p=[];e=[];
	for i in range(len(points)):
	  stat_str = "Testing point %i/%i" % (i,len(points))
	  sys.stdout.write('%s\r' % stat_str)
	  sys.stdout.flush()
	  try:
	    signal.alarm(60)
    	    #print "\tPoint: ", points[i,:]
	    err,comments = estimate_error_2(species_id,points[i,:],tpoints)
	    e.append(err)
	    p.append(points[i,:])
	    #print "\tError: ", err
	    signal.alarm(0)
	    pstr = '%i\t%0.2f\t%0.2f\t%2.2f\t%i\t%2.2f\t%2.2f\t%0.2f\t%i\t%0.2f\t%f\t%s\n' % \
	      (points[i,0], points[i,1], points[i,2], points[i,3], points[i,4], points[i,5], \
		points[i,6], points[i,7], points[i,8], points[i,9], err, comments)
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

def extend_population(population,weights):
    extended_population = []
    num_samples = len(population); num_levels = len(weights)-1
    for i in range(num_samples):
     w_index = i//(num_samples/num_levels)
     for j in range(weights[w_index]):
       extended_population.append(population[i])
    random.shuffle(extended_population)
    return extended_population

def create_children(parent1,parent2,inheritance):
    child1 = parent1; child2 = parent2;
    for i in range(len(parent1)):
      if inheritance[i]==1:
	child1[i]=parent2[i]
	child2[i]=parent1[i]
    return child1, child2

def population_breeding(population,numchildren):
    next_population = []
    for i in range(int(numpy.floor(len(population)/2))):
      for j in range(numchildren):
	parent1 = population[i]; parent2 = population[len(population)-1-i]
	inheritance = numpy.random.choice([0,1], size=(len(parent1)), p=[0.5,0.5])
	child1, child2 = create_children(parent1,parent2,inheritance)
	if numpy.array_equal(child1,child2):
	    next_population.append(child1)
	    continue
	next_population.append(child1)
	next_population.append(child2)
    return next_population

def mutate_child(child,amp,chance):
    mutation_flag = numpy.random.choice([-1,1], size=(len(child)), p=[0.5,0.5])
    for i in range(len(child)):
	if i==0 or i==4:
	  if child[i]<=1:
	    child[i]=child[i]+numpy.abs(mutation_flag[i])
	  else:
	    child[i]=child[i]+mutation_flag[i]
	else:
	    child[i]=child[i]*(1.0+amp*random.random()*mutation_flag[i])
    return child

def mutate_population(population,amp,chance):
    mutated_population=[]
    for i in range(len(population)):
      mutated_population.append(mutate_child(population[i],amp,chance))
    return mutated_population

def check_population(population,ranges):
    # Remove duplicates
    distinct = []
    for i in population:
      if not any(numpy.array_equal(i,j) for j in distinct):
	  distinct.append(i)
    population = distinct;
    clean_population=[]
    for i in range(len(population)):
      locs_low = population[i]<ranges[:,0]
      locs_high = population[i]>ranges[:,1]
      if numpy.any(locs_low):
	population[i][locs_low] = ranges[locs_low,0]
      if numpy.any(locs_high):
	population[i][locs_high] = ranges[locs_high,0]
      clean_population.append(population[i])
    return clean_population


def main():
  return

if __name__ == "__main__":
    main()
