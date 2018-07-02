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
mtg = axialtree2mtg(axialtree, scale, scene)
print list(mtg.property_names())
properties = [(p, 'REAL') for p in mtg.property_names() if p not in ['edge_type', 'index', 'label']]
print properties
mtg_lines = write_mtg(mtg, properties)
print mtg_lines
