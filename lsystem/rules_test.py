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


parameter_dict={'age':' ','trunk_pitch_angle':' ','trunk_roll_angle':' ','trunk_height':' ',
                'no_first_ord_branches':' ','no_second_ord_branches':' ',
                'branching_pitch_angle':' ','branching_roll_angle':' ',
                'diameter_growth_rate':' ','annual_no_new_nodes':' ','avg_internode_length':' '}
parameter_dict['age']=10
parameter_dict['trunk_pitch_angle'] = 5.0            #pitch down wrt turtle's left, in degrees
parameter_dict['trunk_roll_angle'] = 0.0             #roll left wrt turtle's head, in degrees
parameter_dict['trunk_height'] = 3.0            #unit: m, trunk's actual length (regardless of orientation wrt ground) - when the trunk reach this height, it will signal the tree to branch out for the first time
parameter_dict['no_first_ord_branches'] = 3   #number of first order branches
parameter_dict['no_second_ord_branches'] = 5  #number of second order branches (assume uniformity)
parameter_dict['branching_pitch_angle'] = 45.0
parameter_dict['branching_roll_angle'] = 30.0
parameter_dict['diameter_growth_rate'] = 0.1  #m/year? relative growth rate?
parameter_dict['annual_no_new_nodes'] = 30.0  #number of new buds per year
parameter_dict['avg_internode_length'] = 0.03 #unit: m

l = lpy.Lsystem(input_file, parameter_dict)
l.animate()

raw_input() #to prevent window closing at end of program
