
from openalea.mtg import *
from openalea.mtg.io import *
from vplants.plantgl.all import *
import openalea.lpy

g = MTG('C:/Users/Desktop/Desktop/noylum2.mtg')
drf = 'C:/Users/Desktop/Desktop/walnut.drf'
dressing_data = dresser.dressing_data_from_file(drf)
pf = plantframe.PlantFrame(g, TopDiameter='TopDia', DressingData=dressing_data)
pf.plot(gc=True)
Viewer.frameGL.saveImage('C:/Users/Desktop/Desktop/plantframe_example_plantframe-output.png', 'png')

l = openalea.lpy.Lsystem()
axialtree = mtg2axialtree(g)
l.plot(axialtree)
print axialtree
Viewer.frameGL.saveImage('C:/Users/Desktop/Desktop/plantframe_example_lsystem-output.png', 'png')

scale = {'P':1,'A':1,'S':1,'U':1}
scene = l.generateScene(axialtree)
print 'Scene:'
print scene
mtg = axialtree2mtg(axialtree, scale, scene)
print 'Property names:'
print list(mtg.property_names())
properties = [(p, 'REAL') for p in mtg.property_names()]
print 'Properties: '
print properties

import os
current_path = os.path.dirname(os.path.abspath(__file__))
output_mtg = os.path.join(current_path, 'plantframe-example_output.mtg')
mtg_lines = write_mtg(mtg, properties)
f = open(output_mtg, 'w')
f.write(mtg_lines)
f.close()
