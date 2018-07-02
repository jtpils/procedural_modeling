import os
current_path = os.path.dirname(os.path.abspath(__file__))
#input_file = os.path.join(current_path, 'parameters.lpy')
input_file = os.path.join(current_path, 'example.lpy')
output_lpy = os.path.join(current_path, 'lpy_output.png')
output_mtg_topdia = os.path.join(current_path, 'mtg_TopDia_output.png')
output_mtg_diam = os.path.join(current_path, 'mtg_diam_output.png')
output_mtg = os.path.join(current_path, 'mtg_output.png')
mtg_file = os.path.join(current_path, 'lstring_output.mtg')
mtg2d_file = os.path.join(current_path, 'mtg2d.png')
#mtg3d_file = os.path.join(current_path, 'mtg3d.png')

import time
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import numpy
import openalea.lpy as lpy
from openalea.mtg.plantframe import *
from openalea.mtg.mtg import *
from openalea.plantgl.all import *
from openalea.mtg.io import *
from openalea.mtg.aml import *
from openalea.mtg.util import *

l = lpy.Lsystem(input_file)
#l = lpy.Lsystem(input_file, parameter_dict)
axialtree = l.animate()

l.plot(axialtree)
Viewer.frameGL.saveImage(output_lpy, 'png')

scale = {'F':1,'X':1}
#scale = {'A':1,'B':1, 'L':1, 'I':1}
#scene = l.generateScene(axialtree)
scene = l.sceneInterpretation(axialtree)
#scene = l.Tree2Scene(axialtree)
#parameters = {'A':['t', 'o'], 'B':['t', 'o', 'idx'], 'L':['t', 'n'], 'I':['s', 'r']}

#mtg = lpy2mtg(axialtree, l, scene)
#mtg_lines = lpy2mtg(mtg, axialtree, l)
#mtg_lines = write_mtg(mtg)

mtg = axialtree2mtg(axialtree, scale, scene)
#mtg = axialtree2mtg(axialtree, scale, scene, parameter_dict)
#mtg = axialtree2mtg(axialtree, scale, scene, parameters)
#mtg = read_lsystem_string(str(axialtree), scale)
#plot2d(mtg, mtg2d_file, scale)
#plot3d(mtg)

#print "axialtree"
#print axialtree

print list(mtg.property_names())
#properties = [(p, 'REAL') for p in mtg.property_names() if p not in ['edge_type', 'index', 'label']]
properties = [(p, 'REAL') for p in mtg.property_names()]
f = open(mtg_file, 'w')
f.write(write_mtg(mtg, properties))
#f.write(write_mtg(g=mtg, properties=properties, class_at_scale=scale))

#use dict to access mtg structure....
f.close()
