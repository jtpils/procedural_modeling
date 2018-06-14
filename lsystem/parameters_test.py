import os


import numpy as np
#from openalea.pylsystems import *
#os.chdir("D:\Dropbox\Projects\VS1\lpymaster\share\examples")   - use relative path instead
current_path = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.join(current_path, 'parameters.lpy')
import openalea.lpy as lpy
from openalea.lpy import AxialTree
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import time
from openalea.plantgl.all import *
from openalea.mtg.io import lpy2mtg, mtg2lpy, axialtree2mtg, mtg2axialtree
from openalea.mtg.aml import *


parameter_dict={'roll_angle_inc':' ', 'pitch_angle':' ', 'age':' ', 'nb1':' ', 'nb2':' '}
parameter_dict['roll_angle_inc']=60
parameter_dict['pitch_angle']=50
parameter_dict['age']=5
parameter_dict['nb1']=8
parameter_dict['nb2']=np.array([1, 1, 1, 1, 1, 3, 2, 4])

l = lpy.Lsystem(input_file, parameter_dict)
l.animate()

raw_input() #to prevent window closing at end of program
