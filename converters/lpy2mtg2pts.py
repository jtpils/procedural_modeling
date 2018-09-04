import os
import time
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import numpy
import openalea.lpy as lpy
from openalea.mtg.plantframe import *
from openalea.mtg.mtg import *
from openalea.plantgl.all import *
#from openalea.mtg.io import lpy2mtg, mtg2lpy, axialtree2mtg, mtg2axialtree, read_lsystem_string
from openalea.mtg.io import *
from openalea.mtg.aml import *
from openalea.mtg.util import *
from openalea.plantgl.scenegraph._pglsg import *

def lsystem_run(age=10,
                trunk_pitch_angle=5.0, trunk_roll_angle=0.0, trunk_height=3.0,
                no_first_ord_branches=3, no_second_ord_branches=5,
                branching_pitch_angle=45.0, branching_roll_angle=30.0,
                diameter_growth_rate=0.1, annual_no_new_nodes=30.0, avg_internode_length=0.03):
    'Pass known parameter values into L-system rules to produce Lstring output'

    current_path = os.path.dirname(os.path.abspath(__file__))
    #input_file = os.path.join(current_path, 'parameters.lpy')
    input_file = os.path.join(current_path, '../lsystem/rules.lpy')
    input_file = os.path.abspath(os.path.realpath(input_file))
    #input_file = os.path.join(current_path, 'example.lpy')
    output_lpy = os.path.join(current_path, 'lpy_output.png')
    output_mtg_topdia = os.path.join(current_path, 'mtg_TopDia_output.png')
    output_mtg_diam = os.path.join(current_path, 'mtg_diam_output.png')
    output_mtg = os.path.join(current_path, 'mtg_output.png')
    mtg_file = os.path.join(current_path, 'lstring_output.mtg')
    mtg2d_file = os.path.join(current_path, 'mtg2d.png')
    #mtg3d_file = os.path.join(current_path, 'mtg3d.png')
    test_mtg_file = os.path.join(current_path, 'myMtg.mtg')

    parameter_dict={'age':' ','trunk_pitch_angle':' ','trunk_roll_angle':' ','trunk_height':' ',
                    'no_first_ord_branches':' ','no_second_ord_branches':' ',
                    'branching_pitch_angle':' ','branching_roll_angle':' ',
                    'diameter_growth_rate':' ','annual_no_new_nodes':' ','avg_internode_length':' '}
    parameter_dict['age']=age
    parameter_dict['trunk_pitch_angle'] = trunk_pitch_angle            #pitch down wrt turtle's left, in degrees
    parameter_dict['trunk_roll_angle'] = trunk_roll_angle             #roll left wrt turtle's head, in degrees
    parameter_dict['trunk_height'] = trunk_height            #unit: m, trunk's actual length (regardless of orientation wrt ground) - when the trunk reach this height, it will signal the tree to branch out for the first time
    parameter_dict['no_first_ord_branches'] = no_first_ord_branches   #number of first order branches
    parameter_dict['no_second_ord_branches'] = no_second_ord_branches  #number of second order branches (assume uniformity)
    parameter_dict['branching_pitch_angle'] = branching_pitch_angle
    parameter_dict['branching_roll_angle'] = branching_roll_angle
    parameter_dict['diameter_growth_rate'] = diameter_growth_rate  #m/year? relative growth rate?
    parameter_dict['annual_no_new_nodes'] = annual_no_new_nodes  #number of new buds per year
    parameter_dict['avg_internode_length'] = avg_internode_length #unit: m

    #l = lpy.Lsystem(input_file)
    l = lpy.Lsystem(input_file, parameter_dict)
    axialtree = l.animate()
    #axialtree = l.iterate()

    l.plot(axialtree)
    Viewer.frameGL.saveImage(output_lpy, 'png')

    #print axialtree
    #for num, element in enumerate(axialtree):
    #    print num, element

    #scale = {'F':1,'X':1}
    #scale = {'A':1, 'B':1, 'L':1, 'I':1}
    scale = {'I':1}
    #scene = l.generateScene(axialtree)
    scene = l.sceneInterpretation(axialtree)
    #scene = l.Tree2Scene(axialtree)
    #print 'Homomorphism:',
    #print l.homomorphism(axialtree)

    #scene = l.Tree2Scene(axialtree)
    '''
    print 'Scene:'
    print scene
    print 'Components in scene:'
    for shape in scene:
        print shape.geometry.name, shape.geometry.getPglReferenceCount
        #print isinstance(shape, scenegraph._pglsg.Cylinder)
        if type(shape.geometry) is Cylinder:
            print shape.geometry.radius, shape.geometry.height
        elif type(shape.geometry) is Translated:
            print shape.geometry.translation
        elif type(shape.geometry) is Frustum:
            print shape.geometry.radius, shape.geometry.height, shape.geometry.taper
        elif type(shape.geometry) is TriangleSet:
            print 'TriangleSet geometry'
        elif type(shape.geometry) is Oriented:
            print 'Oriented geometry'
        else:
            print 'Unknown shape geometry'
    '''

    #mtg = lpy2mtg(axialtree, l, scene)
    #mtg_lines = lpy2mtg(mtg, axialtree, l)
    #mtg_lines = write_mtg(mtg)

    mtg = axialtree2mtg(axialtree, scale, scene)
    #mtg = axialtree2mtg(axialtree, scale, scene, parameter_dict)
    #parameters = {'A':['age', 'order'], 'B':['age', 'order', 'index'], 'L':['age', 'phyllotactic'], 'I':['order', 'length', 'radius']}
    #parameters = {'A':['t', 'o'], 'B':['t', 'o', 'idx'], 'L':['t', 'n'], 'I':['o', 'a', 'r']}
    #parameters = {'I':['order', 'distanceFromBase', 'length', 'radius']}
    #mtg = axialtree2mtg(axialtree, scale, scene, parameters)

    #extract 3D positions from mtg or scene here! put them into output_point_list
    #todo....

    #mtg = read_lsystem_string(str(axialtree), scale)
    #plot2d(mtg, mtg2d_file, scale)
    #plot3d(mtg)

    #print "axialtree"
    #print axialtree

    #print mtg

    #properties = [(p, 'REAL') for p in mtg.property_names() if p in ['XX', 'YY', 'ZZ']]
    #print properties
    #print write_mtg(mtg, properties)

    #f = open(test_mtg_file)
    #txt = f.read()
    #mtg_test = read_mtg(txt)
    '''
    dressing_data = DressingData(DiameterUnit=1)


    pf1 = PlantFrame(mtg,
                    TopDiameter='TopDia',
                    DressingData = dressing_data)
    #print "TopDia"
    #print pf1.points
    pf1.plot(gc=True)
    Viewer.frameGL.saveImage(output_mtg_topdia, 'png')

    pf2 = PlantFrame(mtg,
                    TopDiameter='diam')
    #print "diam"
    #print pf2.points
    pf2.plot()
    Viewer.frameGL.saveImage(output_mtg_diam, 'png')

    topdia = lambda x:  mtg.property('TopDia').get(x)
    pf3 = PlantFrame(mtg, TopDiameter=topdia, DressingData = dressing_data)
    #axes = pf3._compute_axes(mtg, 3, pf3.points, pf3.origin)
    #axes = pf3._compute_axes()
    #diameters = pf3.algo_diameter()
    #pf3.run()          #called compute_axes & algo_diameter
    #scene = pf3.build_scene(pf3.g, pf3.origin, axes, pf3.points, diameters, 10000) #obsolete function
    #Viewer.display(scene)      #called in plot()
    print "pf3"
    print pf3.points
    pf3.plot(gc=True)
    Viewer.frameGL.saveImage(output_mtg, 'png')
    '''

    #properties = [(p, 'REAL') for p in mtg.property_names() if p not in ['edge_type', 'index', 'label', '_axial_id', 'geometry']]
    properties = [(p, 'REAL') for p in mtg.property_names() if p in ['XX', 'YY', 'ZZ', 'radius']]
    #properties = [(p, 'REAL') for p in mtg.property_names()]
    #properties = [(p, 'REAL') for p in mtg_test.property_names() if p in ['geometry']]
    #print properties
    f = open(mtg_file, 'w')
    #f.write(mtg_lines)
    #f.write(write_mtg(g=mtg, properties=properties, class_at_scale=scale))
    f.write(write_mtg(mtg, properties))

    f.close()

    lstring_output = ''
    return lstring_output


if __name__ == "__main__":
    lsystem_run()
