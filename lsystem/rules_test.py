from __future__ import division
import os
import numpy as np
current_path = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.join(current_path, 'rules.lpy')
import openalea.lpy as lpy
from openalea.lpy import AxialTree
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import time
from openalea.plantgl.all import *
from openalea.mtg.io import lpy2mtg, mtg2lpy, axialtree2mtg, mtg2axialtree
from openalea.mtg.aml import *


parameter_dict={'age':' '}
parameter_dict['age']=1

l = lpy.Lsystem(input_file, parameter_dict)
l.animate()

raw_input() #to prevent window closing at end of program
