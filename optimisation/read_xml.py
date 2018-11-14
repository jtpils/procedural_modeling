import xml.etree.ElementTree as ET
import numpy

def rxml_growthspace(fname):
	tree = ET.parse(fname)
	root = tree.getroot()
	for gs in root.findall('GROWTHSPC'):
		for bgs in gs.findall('B_VOXELS'):
			b_voxsz = float(bgs.find('VOX_SZ').text)
			b_minpt = float(bgs.find('MIN_PT').text)
			b_coords = bgs.find('COORD').text
			b_coords = b_coords.replace('\t','').replace('\n',',').split(',')
			b_coords = numpy.asarray(filter(None,b_coords),dtype=float)
			b_coords = numpy.reshape(b_coords,(-1,3))
		for cgs in gs.findall('C_VOXELS'):
			c_voxsz = float(cgs.find('VOX_SZ').text)
			c_minpt = float(cgs.find('MIN_PT').text)
			c_coords = cgs.find('COORD').text
			c_coords = c_coords.replace('\t','').replace('\n',',').split(',')
			c_coords = numpy.asarray(filter(None,c_coords),dtype=float)
			c_coords = numpy.reshape(c_coords,(-1,3))
	return b_voxsz, b_voxsz*b_coords, c_voxsz, c_voxsz*c_coords

def rxml_treeparams(fname):
	tree = ET.parse(fname)
	root = tree.getroot()
	for identity in root.findall('IDENTITY'):
		cname = identity.find('COMMON').text
		species = identity.find('BOTANICAL').text
		location = identity.find('BASE_COORD').text
		location = numpy.asarray(location.split(","),dtype=float)
	for vals in root.findall('SIZE'):
		height = float(vals.find('HEIGHT').text)
	trunk_height = 0.0
	for branch in root.findall('BRANCHING'):
		max_order = float(branch.find('ORDER').text)
		for seg in branch.findall('SEG'):
		  segment = seg.find('S_CORD').text
		  segment = numpy.asarray(segment.split(','),dtype=float)
		  trunk_height = trunk_height + segment[2]
	return cname, species, location, height, trunk_height

def main():
	xml_fname = '../../growth-space/example_xml/Tree1_Parameters.xml'
	common_name, species, location, height, trunk_height = rxml_treeparams(xml_fname)
	b_voxsz, b_coords, c_voxsz, c_coords = rxml_growthspace(xml_fname)


	print "#-------------------------------------------------#"
	print "#   Common name:\t%s" % common_name
	print "#   Botanic name:\t%s" % species
	print "#   Location"
	print "#   \tx:\t\t%f\n#\ty:\t\t%f\n#\tz:\t\t%f" % (location[0],location[1],location[2])
	print "#   Tree height:\t%f" % height
	print "#   Trunk height:\t%f" % trunk_height
	print "#   Voxel sizes"
	print "#   \tb:\t\t%f\n#\tc:\t\t%f" % (b_voxsz, c_voxsz)
	print "#   Number of voxels"
	print "#   \tb:\t\t%i\n#\tc:\t\t%i" % (len(b_coords),len(c_coords))
	print "#   Bounding box:\tbranch"
	print "#   \tx:\t\t(%f) - (%f)" % (numpy.min(b_coords[:,0]),numpy.max(b_coords[:,0]))
	print "#   \ty:\t\t(%f) - (%f)" % (numpy.min(b_coords[:,1]),numpy.max(b_coords[:,1]))
	print "#   \tz:\t\t(%f) - (%f)" % (numpy.min(b_coords[:,2]),numpy.max(b_coords[:,2]))
	print "#   Bounding box:\tcrown"
	print "#   \tx:\t\t(%f) - (%f)" % (numpy.min(c_coords[:,0]),numpy.max(c_coords[:,0]))
	print "#   \ty:\t\t(%f) - (%f)" % (numpy.min(c_coords[:,1]),numpy.max(c_coords[:,1]))
	print "#   \tz:\t\t(%f) - (%f)" % (numpy.min(c_coords[:,2]),numpy.max(c_coords[:,2]))
	print "#-------------------------------------------------#"

if __name__ == "__main__":
	main()
