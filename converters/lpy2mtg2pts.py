#Note: to copy file 'io-move.py' in this folder to the equivalent location of 'C:\Python27\Lib\site-packages\OpenAlea.Mtg-1.x.0-py2.7.egg\openalea\mtg' (rename to io.py and overwrite existing io.py to modify mtg writing)

import os
from enum import Enum
from datetime import datetime
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import openalea.lpy as lpy
from openalea.mtg.plantframe import *
from openalea.mtg.mtg import *
from openalea.plantgl.all import *
from openalea.mtg.io import *
from openalea.mtg.aml import *
from openalea.mtg.util import *
from openalea.plantgl.scenegraph._pglsg import *
import cStringIO
import optimisation.read_xml as rxml

class Species(Enum):
    Unspecified = 0
    AA = 1  #Archontophoenix alexandrae (palm)
    SS = 2  #Samanea saman (raintree)
    PP = 3  #Peltophorum pterocarpum (yellow flame)
    HO = 4  #Hopea odorata
    SMa = 5  #Swietenia macrophylla (mahogany)
    KS = 6  #Khaya senegalensis
    SG = 7  #Syzygium grande
    TR = 8  #Tabebuia rosea
    SMy = 9  #Syzygium myrtifolium
    SP = 10  #Sterculia parviflora


def lsystem_run(species=Species.Unspecified, age=10,
                trunk_pitch_angle=5.0, trunk_roll_angle=0.0, trunk_height=3.0,
                no_first_ord_branches=3,
                branching_pitch_angle=45.0, branching_roll_angle=30.0,
                diameter_growth_rate=0.1, annual_no_new_nodes=30.0, avg_internode_length=0.03):
    'Pass known parameter values into L-system rules to produce Lstring output'
    flag_printString = False
    flag_animate = True
    flag_plot = False
    flag_writeToFile = False

    current_path = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_path, '../lsystem/rules.lpy')
    input_file = os.path.abspath(os.path.realpath(input_file))
    output_lpy = os.path.join(current_path, 'lpy_output.png')
    mtg_file = os.path.join(current_path, './mtg/'+str(datetime.now().strftime('%Y%m%d%H%M%S'))+'_lstring_'+str(species)+'-age'+str(age)+'-trunkpitch'+str(trunk_pitch_angle)+',roll'+str(trunk_roll_angle)+',height'+str(trunk_height)+'-1ordbranches'+str(no_first_ord_branches)+'-branchpitch'+str(branching_pitch_angle)+',roll'+str(branching_roll_angle)+'-diamgrowth'+str(diameter_growth_rate)+'-ann_nodes'+str(annual_no_new_nodes)+'-intlen'+str(avg_internode_length)+'.mtg')

    parameter_dict = {'species':' ','age':' ','trunk_pitch_angle':' ','trunk_roll_angle':' ','trunk_height':' ',
                      'no_first_ord_branches':' ',
                      'branching_pitch_angle':' ','branching_roll_angle':' ',
                      'diameter_growth_rate':' ','annual_no_new_nodes':' ','avg_internode_length':' '}
    parameter_dict['species'] = species
    parameter_dict['age'] = age
    parameter_dict['trunk_pitch_angle'] = trunk_pitch_angle            #pitch down wrt turtle's left, in degrees
    parameter_dict['trunk_roll_angle'] = trunk_roll_angle             #roll left wrt turtle's head, in degrees
    parameter_dict['trunk_height'] = trunk_height            #unit: m, trunk's actual length (regardless of orientation wrt ground) - when the trunk reach this height, it will signal the tree to branch out for the first time
    parameter_dict['no_first_ord_branches'] = no_first_ord_branches   #number of first order branches
    parameter_dict['branching_pitch_angle'] = branching_pitch_angle
    parameter_dict['branching_roll_angle'] = branching_roll_angle
    parameter_dict['diameter_growth_rate'] = diameter_growth_rate  #r of diameter growth rate RGR (relative growth rate) equation model
    parameter_dict['annual_no_new_nodes'] = annual_no_new_nodes  #number of new buds per year
    parameter_dict['avg_internode_length'] = avg_internode_length #unit: m

    lsys = lpy.Lsystem(input_file, parameter_dict)
    #lsys.useGroup(0)       #use rules for determinate growth pattern

    if flag_animate is True:
        axialtree = lsys.animate()
    else:
        axialtree = lsys.iterate()  #lstring output muted in rules.lpy's function EndEach

    if flag_printString is True:
        print 'L-string output:'
        print axialtree

    if flag_plot is True:
        lsys.plot(axialtree)
        Viewer.frameGL.saveImage(output_lpy, 'png')

    scale = {'I':1}
    scene = lsys.sceneInterpretation(axialtree)
    mtg = axialtree2mtg(axialtree, scale, scene)

    properties = [(p, 'REAL') for p in mtg.property_names() if p in ['XX', 'YY', 'ZZ', 'radius']]

    mtg_string_output = write_mtg(mtg, properties)

    if flag_writeToFile is True:
        f = open(mtg_file, 'w')
        f.write(mtg_string_output)
        f.close()

    return mtg_string_output


def mtg2obj(mtg_string):
    output = cStringIO.StringIO()
    output.write('#This OBJ file was generated from MTG file.\n')

    buffer = cStringIO.StringIO(mtg_string)
    for line in buffer:
        if line.startswith(('/', '+', '^', '\t')):
            pieces = line.split('\t')
            print >>output, 'v', pieces[len(pieces)-3], pieces[len(pieces)-4], pieces[len(pieces)-1],

    obj_string = output.getvalue()
    output.close()

    return obj_string


def test_mtg2obj():
    current_path = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_path, './mtg/test.mtg')
    input_file = os.path.abspath(os.path.realpath(input_file))
    output_file = os.path.join(current_path, './obj/'+str(datetime.now().strftime('%Y%m%d%H%M%S'))+'_test_mtg.obj')

    with open(input_file, 'r') as in_file:
        mtg_string = in_file.read()

    obj_string = mtg2obj(mtg_string)

    f = open(output_file, 'w')
    f.write(obj_string)
    f.close()


def xmlgrowthspace2obj(xml_filename):
    output = cStringIO.StringIO()
    output.write('#This OBJ file was generated from XML growth space file.\n')

    voxel_size, coords = rxml.rxml_growthspace(xml_filename)

    for point in coords:
        print >>output, 'v', voxel_size*point[0], voxel_size*point[1], voxel_size*point[2]

    obj_string = output.getvalue()
    output.close()

    return obj_string


def test_xmlgrowthspace2obj():
    current_path = os.path.dirname(os.path.abspath(__file__))
    input_file = os.path.join(current_path, './xml/test.xml')
    input_file = os.path.abspath(os.path.realpath(input_file))
    output_file = os.path.join(current_path, './obj/'+str(datetime.now().strftime('%Y%m%d%H%M%S'))+'_test_xmlgrowthspace.obj')

    #with open(input_file, 'r') as in_file:
    #        xml_string = in_file.read()

    obj_string = xmlgrowthspace2obj(input_file)

    f = open(output_file, 'w')
    f.write(obj_string)
    f.close()


if __name__ == "__main__":
    lsystem_run(species=Species.SS, age=25,
                trunk_pitch_angle=-0.39, trunk_roll_angle=-0.37, trunk_height=7.67,
                no_first_ord_branches=0,
                branching_pitch_angle=0.0, branching_roll_angle=0.0,
                diameter_growth_rate=0.6125, annual_no_new_nodes=7, avg_internode_length=0.1537)
    print 'Program finished'
    raw_input()

    #test_mtg2obj()
    #test_xmlgrowthspace2obj()
