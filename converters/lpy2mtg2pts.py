import os
current_path = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.join(current_path, 'example.lpy')
output_file = os.path.join(current_path, 'example_output.png')
mtg_file = os.path.join(current_path, 'lstring_output.mtg')

import openalea.lpy as lpy
from openalea.lpy import AxialTree
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import time
from openalea.plantgl.all import *
from openalea.mtg.io import lpy2mtg, mtg2lpy, axialtree2mtg, mtg2axialtree, read_lsystem_string
from openalea.mtg.aml import *

l = lpy.Lsystem(input_file)

axialtree = l.iterate()

l.plot(axialtree)

Viewer.frameGL.saveImage(output_file, 'png')

scale = {'F':1,'X':1}
#mtg = lpy2mtg(axialtree, l)
#mtg_lines = lpy2mtg(mtg, axialtree, l)
#mtg_lines = write_mtg(mtg)

mtg = read_lsystem_string(str(axialtree), scale)

f = open(mtg_file, 'w')
#f.write(mtg_lines)
f.write(str(mtg))

#use dict to access mtg structure....

f.close()