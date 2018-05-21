import os


import numpy as np
#from openalea.pylsystems import *
os.chdir("D:\Dropbox\Projects\VS1\lpymaster\share\examples")
import openalea.lpy as lpy
from openalea.lpy import AxialTree
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import time
from openalea.plantgl.all import *
from openalea.mtg.io import lpy2mtg, mtg2lpy, axialtree2mtg, mtg2axialtree
from openalea.mtg.aml import *


parameter_dict={'roll_angle_inc':' ', 'pitch_angle':' ', 'age':' ', 'nb1':' ', 'nb2':' '}
parameter_dict['roll_angle_inc']=120
parameter_dict['pitch_angle']=30
parameter_dict['age']=4
parameter_dict['nb1']=4
parameter_dict['nb2']=np.array([2, 1, 3, 4])

l = lpy.Lsystem("parameters.lpy", parameter_dict)
l.animate()

