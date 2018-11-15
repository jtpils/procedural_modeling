#from __future__ import division
import os
current_path = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.join(current_path, 'rules.lpy')
import openalea.lpy as lpy
#from openalea.lpy import AxialTree
#from PyQt4.QtCore import *
#from PyQt4.QtGui import *
#import time
#from openalea.plantgl.all import *
#from openalea.mtg.io import lpy2mtg, mtg2lpy, axialtree2mtg, mtg2axialtree
#from openalea.mtg.aml import *
from enum import Enum

class Species(Enum):
    Undefined = 0
    AA = 1  #Archontophoenix alexandrae (palm)
    SS = 2  #Samanea saman (raintree)
    PP = 3  #Peltophorum pterocarpum (yellow flame)
    HO = 4  #Hopea odorata
    SMa = 5  #Swietenia macrophylla (mahogany)
    KS = 6  #Khaya senegalensis
    SG = 7  #Syzygium grande
    TR = 8  #Tabebuia rosea
    SMy = 9  #Syzygium myrtifolium
    SP = 10  #Sterculia parviflora

parameter_dict={'species':' ','age':' ','trunk_pitch_angle':' ','trunk_roll_angle':' ','trunk_height':' ',
                'no_first_ord_branches':' ','no_second_ord_branches':' ',
                'branching_pitch_angle':' ','branching_roll_angle':' ',
                'diameter_growth_rate':' ','annual_no_new_nodes':' ','avg_internode_length':' '}
parameter_dict['species']=Species.SS
parameter_dict['age']=30
parameter_dict['trunk_pitch_angle'] = 0.0            #pitch down wrt turtle's left, in degrees
parameter_dict['trunk_roll_angle'] = 0.0             #roll left wrt turtle's head, in degrees
parameter_dict['trunk_height'] = 1.0            #unit: m, trunk's actual length (regardless of orientation wrt ground) - when the trunk reach this height, it will signal the tree to branch out for the first time
parameter_dict['no_first_ord_branches'] = 1   #number of first order branches
parameter_dict['branching_pitch_angle'] = 45.0
parameter_dict['branching_roll_angle'] = 120.0
parameter_dict['diameter_growth_rate'] = 0.01588  #m/year? relative growth rate?
parameter_dict['annual_no_new_nodes'] = 24.0  #number of new buds per year
parameter_dict['avg_internode_length'] = 0.02137 #unit: m

l = lpy.Lsystem(input_file, parameter_dict)
l.animate()

raw_input() #to prevent window closing at end of program
