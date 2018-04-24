import os
current_path = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.join(current_path, 'example.lpy')
output_file = os.path.join(current_path, 'example_output.png')
mtg_file = os.path.join(current_path, 'lstring_output.mtg')
mtg2d_file = os.path.join(current_path, 'mtg2d.png')
#mtg3d_file = os.path.join(current_path, 'mtg3d.png')

import openalea.lpy as lpy
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import time
from openalea.plantgl.all import *
#from openalea.mtg.io import lpy2mtg, mtg2lpy, axialtree2mtg, mtg2axialtree, read_lsystem_string
from openalea.mtg.io import *
from openalea.mtg.aml import *
from openalea.mtg.util import *

l = lpy.Lsystem(input_file)

axialtree = l.iterate()

l.plot(axialtree)

Viewer.frameGL.saveImage(output_file, 'png')

scale = {'F':1,'X':1}
scene = l.sceneInterpretation(axialtree)

#mtg = lpy2mtg(axialtree, l)
#mtg_lines = lpy2mtg(mtg, axialtree, l)
#mtg_lines = write_mtg(mtg)

axialtree2mtg(axialtree, scale, scene)
mtg = read_lsystem_string(str(axialtree), scale)
plot2d(mtg, mtg2d_file)
plot3d(mtg)

f = open(mtg_file, 'w')
#f.write(mtg_lines)
f.write(write_mtg(g=mtg, class_at_scale=scale))



#use dict to access mtg structure....

f.close()
