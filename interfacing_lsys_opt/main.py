# Created by Like on 7 Mar 2018

import numpy
from pyquaternion import Quaternion
import os
import sys
sys.path.append('../converters/')
#from converters.lpy2mtg2pts import lsystem_run
from lpy2mtg2pts import lsystem_run

def estimate_age(species, size):
    'Estimate tree age based on scanned size of tree of certain species'
    #todo - estimate based on species, lookup table for tree size/height/dimension vs age range
    age = 4
    return age

def parse_lstring(lstring):
    #todo

    return output_point_list

def generate_lsystem_tree_points(params):
    'Wrapper function to map params array to individual variables'
    'Note: all parameters are optional with default values'

    return ngenerate_lsystem_tree_points(   age=int(params[0]),
                                            trunk_pitch_angle=params[1],
                                            trunk_roll_angle=params[2],
                                            trunk_height=params[3],
                                            no_first_ord_branches=int(params[4]),
                                            no_second_ord_branches=int(params[5]),
                                            branching_pitch_angle=params[6],
                                            branching_roll_angle=params[7],
                                            diameter_growth_rate=params[8],
                                            annual_no_new_nodes=params[9],
                                            avg_internode_length=params[10] )

def ngenerate_lsystem_tree_points(  age=1,
                                    trunk_pitch_angle=5.0,
                                    trunk_roll_angle=0.0,
                                    trunk_height=3.0,
                                    no_first_ord_branches=3,
                                    no_second_ord_branches=5,
                                    branching_pitch_angle=45.0,
                                    branching_roll_angle=30.0,
                                    diameter_growth_rate=0.1,
                                    annual_no_new_nodes=30.0,
                                    avg_internode_length=0.03 ):

    'Parse L-System string output into a list of points of cylinder base centers and their radii forming the trunk-branch skeleton representation of a tree'

    output_point_list = lsystem_run(age, trunk_pitch_angle, trunk_roll_angle, trunk_height,
                                    no_first_ord_branches, no_second_ord_branches,
                                    branching_pitch_angle, branching_roll_angle,
                                    diameter_growth_rate, annual_no_new_nodes, avg_internode_length)
    # return parse_lstring(lstring)
    #parse l-system output string into 3d points matching Xfrog object size

    return output_point_list

def estimate_error(output_point_list, growth_space):

    return cost;

def optimise():
    age = estimate_age()
    cost = 0
    threshold = 100

    #parameters (with their PDFs for species X) to optimise
    no_1st_ord_branches = 3 #valid range [2,4]
    no_2nd_ord_branches = 2 #valid range [2,4]
    branching_angle_roll = 0.0 #valid range [0,360]
    branching_angle_pitch = 20.0 #valid range [10,90]

    if (cost>threshold):
        #backtrack
        pass
    else:  #
        #todo - decide on the values of parameters by sampling the PDF
        output_point_list = generate_lsystem_tree(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch)
        cost = estimate_error(output_point_list, growth_space)
    return growth_param_list

def main():
    #output = generate_lsystem_tree_points(numpy.array([10, 5.0, 0.0, 3.0, 3, 5, 45.0, 30.0, 0.1, 30.0, 0.03]))
    output = ngenerate_lsystem_tree_points( age=10,
                                            trunk_pitch_angle=5.0,
                                            trunk_roll_angle=0.0,
                                            trunk_height=3.0,
                                            no_first_ord_branches=3,
                                            no_second_ord_branches=5,
                                            branching_pitch_angle=45.0,
                                            branching_roll_angle=30.0,
                                            diameter_growth_rate=0.1,
                                            annual_no_new_nodes=30.0,
                                            avg_internode_length=0.03 )
    print output
    file = open('lsys_optimisation_output.vtk', 'w')
    file.write("# vtk DataFile Version 3.0\nTree points\nASCII\nDATASET POLYDATA\n\nPOINTS "+str(len(output))+" integer ")
    for item in output:
        file.write("%d %d %d\n" % (item[0], item[1], item[2]))
    file.close()

if __name__ == "__main__":
    main()
