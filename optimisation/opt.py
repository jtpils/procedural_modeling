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

def estimate_error(sp_id,params,tpts,tg):
  # Estimate error from rasterised tings
  # Collect points
  suppress()
  ls_mtg = interf.generate_lsystem_tree_points(sp_id,params)
  allow()
  ls_pts,ls_rad = lgs.mtg_string_gs(ls_mtg)
  #print ls_pts, ls_rad
  # Compute target vars
  t_bbox, t_th, t_ch, t_mu, t_el, t_eu = cgs.cgs_vx(tpts)
  # Compute LS vars
  ls_bbox, ls_th, ls_ch, ls_mu, ls_el, ls_eu = cgs.cgs_ls(ls_pts)
  # Error in horizontal bbox
  E1 = (numpy.sqrt((t_bbox[1]-t_bbox[0])**2+(t_bbox[3]-t_bbox[2])**2)- \
    numpy.sqrt((ls_bbox[1]-ls_bbox[0])**2+(ls_bbox[3]-ls_bbox[2])**2))/ \
    numpy.sqrt((t_bbox[1]-t_bbox[0])**2+(t_bbox[3]-t_bbox[2])**2)
  #comm_e1 = "Error in hbb=%f" % E1

  # Error in tree height
  E2 = numpy.abs(t_bbox[5]-ls_bbox[5])/t_bbox[5]
  comm_e2 = "Tree heights: GS=%4.2f, LS=%4.2f " % (t_bbox[5],ls_bbox[5])

  E3 = 0.25*( \
    numpy.abs(t_el[1]-ls_el[1])/t_el[1] + numpy.abs(t_el[2]-ls_el[2])/t_el[2] + \
    numpy.abs(t_eu[1]-ls_eu[1])/t_eu[1] + numpy.abs(t_eu[2]-ls_eu[2])/t_eu[2])

  E = E1 + E2 + E3
  comm_E = "| E1=%4.2f, E2=%4.2f, E3=%4.2f" % (E1,E2,E3)

  comments = "# " + comm_e2 + comm_E

  return E, comments


def estimate_error_2(sp_id,params,tpts,tg):
  # Estimate error from rasterised tings
  # Collect points
  suppress()
  ls_mtg = interf.generate_lsystem_tree_points(sp_id,params)
  allow()
  ls_pts,ls_rad = lgs.mtg_string_gs(ls_mtg)
  #ls_pts=ls_pts/10.
  #ls_rad=ls_rad/10.
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

  E3 = 0.5*(numpy.abs((lsbbox[0,1]-lsbbox[0,0])-(tbbox[0,1]-tbbox[0,0]))/(tbbox[0,1]-tbbox[0,0]) + \
    numpy.abs((lsbbox[1,1]-lsbbox[1,0])-(tbbox[1,1]-tbbox[1,0]))/(tbbox[1,1]-tbbox[1,0]))

  E4 = numpy.abs(lsbbox[2,1]-tbbox[2,1])/tbbox[2,1]

  E5 = 0.0
  if (lsbbox[2,0]<0):
    E5 = 1.0

  E6 = numpy.abs(ls_rad-tg)/tg



  #E3 = 0.0

  if (numpy.abs((lsbbox[0,1]-lsbbox[0,0])-(tbbox[0,1]-tbbox[0,0]))/(tbbox[0,1]-tbbox[0,0])>0.5) \
    or (numpy.abs((lsbbox[1,1]-lsbbox[1,0])-(tbbox[1,1]-tbbox[1,0]))/(tbbox[1,1]-tbbox[1,0])>0.5):
    #E3 = E3 + 1.0
    comments = comments + "W "
  else:
    comments = comments + "  "

  if numpy.abs(lsbbox[2,1]-tbbox[2,1])/tbbox[2,1]>0.25:
    #E3 = E3 + 1.0
    comments = comments + "H "
  else:
    comments = comments + "  "

  if (lsbbox[2,0]<0):
    #E3 = E3 + 1.0
    comments = comments + "0 "
  else:
    comments = comments + "  "

  #E4=0.0
  E_trunk = numpy.abs(ls_rad-tg)/tg;
  if E_trunk>0.5:
    #E4=E4+1.
    comments = comments + "T "
  else:
    comments = comments + "  "


  app_str = "| LS_H=%f, GS_H=%f, LS_T=%f, GS_T=%f" % (lsbbox[2,1], tbbox[2,1],ls_rad,tg)
  comments = comments + app_str

  # E1-space, E2-full raster, E3-bbox, E4-height, E5-zero, E6-trunk
  E = 0.15*E1 + 0.05*E2 + 0.2*E3 + 0.2*E4 + 0.2*E5 + 0.2*E6

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

def map_error(species_id,points,tpoints,tg,fname):
	p=[];e=[];
	for i in range(len(points)):
	  stat_str = "Testing point %i/%i" % (i,len(points))
	  sys.stdout.write('%s\r' % stat_str)
	  sys.stdout.flush()
	  try:
	    signal.alarm(180)
    	    #print "\tPoint: ", points[i,:]
	    err,comments = estimate_error_2(species_id,points[i,:],tpoints,tg)
	    e.append(err)
	    p.append(points[i,:])
	    #print "\tError: ", err
	    signal.alarm(0)
	    pstr = '%i\t%0.2f\t%0.2f\t%2.2f\t%i\t%2.4f\t%2.2f\t%0.3f\t%i\t%0.4f\t%f\t%s\n' % \
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
