import read_xml as rxml
import numpy

def xml_file_gs(fname):
	b_voxsz, b_coords, c_voxsz, c_coords = rxml.rxml_growthspace(fname)
	gs = numpy.vstack((b_coords, c_coords))
	return gs

def mtg_file_gs(fname):
	mtg_file = fname
	#print "Reading file: %s" % mtg_file
	with open(mtg_file) as f:
		lines = f.readlines()
	i=0
	while(1):
		if lines[i].startswith("ENTITY-CODE"):
			break
		else:
			i=i+1
	# Columns are 1=y, 2=radius, 3=z, 4=z
	pts = numpy.genfromtxt(mtg_file,skip_header=i+1,usecols=(4,1,3))
	return pts

def normalised_voxel_gs(fname,resolution):
	# Function reads growth_space file, shifts x,y coordinates
    # to zero mean, sets min z coordinate to 0.
    gs = numpy.genfromtxt(fname,delimiter=',')
    gs = gs*resolution/100.
    gs[:,0] = gs[:,0] - numpy.mean(gs[:,0])
    gs[:,1] = gs[:,1] - numpy.mean(gs[:,1])
    gs[:,2] = gs[:,2] - numpy.min(gs[:,2])
    return gs
