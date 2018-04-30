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

import openalea.lpy as lpy
from openalea.mtg.plantframe import *
from openalea.mtg.mtg import *
from openalea.plantgl.all import *
#from openalea.mtg.io import lpy2mtg, mtg2lpy, axialtree2mtg, mtg2axialtree, read_lsystem_string
from openalea.mtg.io import *
from openalea.mtg.aml import *
from openalea.mtg.util import *

l = lpy.Lsystem(input_file)

axialtree = l.iterate()
#axialtree = l.iterate(7)

l.plot(axialtree)
Viewer.frameGL.saveImage(output_lpy, 'png')

scale = {'F':1,'X':1}
#scale = {'A':1,'B':1, 'L':1, 'I':1}
scene = l.sceneInterpretation(axialtree)
#scene = l.Tree2Scene(axialtree)

#mtg = lpy2mtg(axialtree, l, scene)
#mtg_lines = lpy2mtg(mtg, axialtree, l)
#mtg_lines = write_mtg(mtg)

mtg = axialtree2mtg(axialtree, scale, scene)
#mtg = read_lsystem_string(str(axialtree), scale, scene)
#mtg = read_lsystem_string(str(axialtree), scale)
plot2d(mtg, mtg2d_file, scale)
#plot3d(mtg)

dressing_data = DressingData(DiameterUnit=5)
pf1 = PlantFrame(mtg,
                TopDiameter='TopDia',
                DressingData = dressing_data)
print "TopDia"
print pf1.points
pf1.plot(gc=True)
Viewer.frameGL.saveImage(output_mtg_topdia, 'png')

pf2 = PlantFrame(mtg,
                TopDiameter='diam')
print "diam"
print pf2.points
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

f = open(mtg_file, 'w')
#f.write(mtg_lines)
f.write(write_mtg(g=mtg, class_at_scale=scale))
#use dict to access mtg structure....
f.close()
