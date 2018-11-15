# Created by Like on 7 Mar 2018

import numpy
from pyquaternion import Quaternion
import os
import sys
from converters.lpy2mtg2pts import lsystem_run, Species
import warnings
warnings.filterwarnings("ignore")


def estimate_age(species, size):
    'Estimate tree age based on scanned size of tree of certain species'
    #todo - estimate based on species, lookup table for tree size/height/dimension vs age range
    age = 4
    return age


def generate_lsystem_tree_points(species_id, params):
    'Wrapper function to map params array to individual variables'
    'Note: all parameters are optional with default values'

    return ngenerate_lsystem_tree_points(   species=Species(species_id).name,
                                            age=int(params[0]),
                                            trunk_pitch_angle=params[1],
                                            trunk_roll_angle=params[2],
                                            trunk_height=params[3],
                                            no_first_ord_branches=int(params[4]),
                                            branching_pitch_angle=params[5],
                                            branching_roll_angle=params[6],
                                            diameter_growth_rate=params[7],
                                            annual_no_new_nodes=params[8],
                                            avg_internode_length=params[9] )


def ngenerate_lsystem_tree_points(  species=Species.Undefined,
                                    age=1,
                                    trunk_pitch_angle=5.0,
                                    trunk_roll_angle=0.0,
                                    trunk_height=3.0,
                                    no_first_ord_branches=3,
                                    branching_pitch_angle=45.0,
                                    branching_roll_angle=30.0,
                                    diameter_growth_rate=0.1,
                                    annual_no_new_nodes=30.0,
                                    avg_internode_length=0.03 ):

    'Parse L-System string output into a list of points of cylinder base centers and their radii forming the trunk-branch skeleton representation of a tree'

    output_point_list = lsystem_run(species, age, trunk_pitch_angle, trunk_roll_angle, trunk_height,
                                    no_first_ord_branches,
                                    branching_pitch_angle, branching_roll_angle,
                                    diameter_growth_rate, annual_no_new_nodes, avg_internode_length)
    # return parse_lstring(lstring)
    #parse l-system output string into 3d points matching Xfrog object size

    return output_point_list


'''
def estimate_error(output_point_list, growth_space):
    cost = 0
    return cost


def optimise():
    species = Species.Undefined
    age = estimate_age()
    cost = 0
    threshold = 100

    #parameters (with their PDFs for species X) to optimise
    no_1st_ord_branches = 3 #valid range [2,4]
    branching_angle_roll = 0.0 #valid range [0,360]
    branching_angle_pitch = 20.0 #valid range [10,90]

    if (cost > threshold):
        #backtrack
        pass
    else:  #
        #todo - decide on the values of parameters by sampling the PDF
        output_point_list = generate_lsystem_tree(species, age, no_1st_ord_branches, branching_angle_roll, branching_angle_pitch)
        cost = estimate_error(output_point_list, growth_space)
    return growth_param_list
'''


def main():
    #output = generate_lsystem_tree_points(Species.Undefined, numpy.array([10, 5.0, 0.0, 3.0, 3, 45.0, 30.0, 0.1, 30.0, 0.03]))
    output = ngenerate_lsystem_tree_points( species=Species.Undefined,
                                            age=20,
                                            trunk_pitch_angle=0.0,
                                            trunk_roll_angle=0.0,
                                            trunk_height=3.0244,
                                            no_first_ord_branches=3,
                                            branching_pitch_angle=45.0,
                                            branching_roll_angle=120.0,
                                            diameter_growth_rate=0.01588,
                                            annual_no_new_nodes=24.0,
                                            avg_internode_length=0.02137 )

    print output
#    file = open('lsys_optimisation_output.vtk', 'w')
#   file.write("# vtk DataFile Version 3.0\nTree points\nASCII\nDATASET POLYDATA\n\nPOINTS "+str(len(output))+" integer ")
#    for item in output:
#        file.write("%d %d %d\n" % (item[0], item[1], item[2]))
#    file.close()


if __name__ == "__main__":
    main()
