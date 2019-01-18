import optimisation.read_xml as rxml
import numpy

def xml_file_gs(fname,split):
	if split==1:
	  b_voxsz, b_coords, c_voxsz, c_coords = rxml.rxml_seg_growthspace(fname)
	  gs = numpy.vstack((b_coords, c_coords))
	if split==0:
	  voxsz, coords = rxml.rxml_growthspace(fname)
	  gs = coords
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
	# Columns are 1=y, 2=x, 3=radius, 4=z
	pts = numpy.genfromtxt(mtg_file,skip_header=i+1,usecols=(2,1,4))
	return pts

def mtg_string_gs(mtg_str):
	split_str = []
	for line in mtg_str.split('\n'):
		line = line.replace("^","");
		line = line.replace("<","");
		line = line.replace(">","");
		line = line.replace("I","");
		line = line.replace("+","");
		line = line.replace("/","");
		line = line.lstrip()
		new_line = line.replace("\t"," ")
		split_str.append(new_line)
	#print numpy.shape(split_str)
	i=0
	while(1):
		if split_str[i].startswith("ENTTY-CODE"):
			break
		else:
			i=i+1
	#print split_str[i-1]
	# Columns are 1=y, 2=radius, 3=z, 4=z
	pts_list = numpy.asarray(split_str[i+1:-1])

	pts = []
	for j in range(len(pts_list)):
		line = pts_list[j].split()
		num_line = [float(x) for x in line]
		#print line, num_line
		pts.append(num_line)
	pts = numpy.asarray(pts)
 	pts = pts[:,[1,0,3]]
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

def main():
	xml_fname = '../converters/raintree.mtg'
	pts = mtg_file_gs(xml_fname)
	print numpy.min(pts[:,0]), numpy.min(pts[:,1]), numpy.min(pts[:,2])
	print numpy.max(pts[:,0]), numpy.max(pts[:,1]), numpy.max(pts[:,2])

if __name__ == "__main__":
	main()
