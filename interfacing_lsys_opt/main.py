# Created by Like on 7 Mar 2018

import numpy
from pyquaternion import Quaternion

def estimate_age():
    age = 4
    return age

def lsystem_run(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch):
    'Pass known parameter values into L-system rules to produce Lstring output'
    #todo - Peng

    lstring_output = ''

    return lstring_output

def generate_lsystem_tree(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch):
    'Translate Lstring output into a list of points of cylinder base centers and their radii forming the trunk-branch representation of a tree'

    current_post = numpy.array([0,0,0])
    if age == 0:
        return current_post

    out = lsystem_run(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch)
    #parse l-system output string
    #todo
    #to match with Xfrog object size

    # temporary hardcoding of tree growth: age 1-2 trunk grows taller by 2 m/yr, age 3 first branch out by 1 m/yr, age 4 second branch out by 0.5 m/yr
    # ========================
    output_point_list = [[0,0,0]]  #ground
    for t in range (1, age+1): #trunk
        trunk_height = [0,0,2*t]
        output_point_list = output_point_list + [trunk_height]

    v = [0,0,1]
    pitch_axis = [1,0,0]
    roll_axis = [0,0,1]
    pitch_angle = branching_angle_pitch/180 * 3.14 #radian
    roll_angle = branching_angle_roll/180 * 3.14 #radian

    if age > 2:
        first_pitched_v = Quaternion(axis=pitch_axis,angle=pitch_angle).rotate(v)
        for t in range (0, no_1st_ord_branches): #1st order branch
            first_rolled_v = Quaternion(axis=roll_axis,angle=roll_angle).rotate(first_pitched_v)
            output_point_list = output_point_list + [[a+b for a,b in zip(trunk_height,first_rolled_v)]]

            if age > 3:
                new_pitch_axis = Quaternion(axis=roll_axis,angle=roll_angle).rotate(pitch_axis)
                new_roll_axis = first_rolled_v
                second_pitched_v = Quaternion(axis=new_pitch_axis,angle=pitch_angle).rotate([0.5*x for x in first_rolled_v])
                for u in range (0, no_2nd_ord_branches): #2nd order branch
                    second_rolled_v = Quaternion(axis=new_roll_axis,angle=roll_angle).rotate(second_pitched_v)
                    output_point_list = output_point_list + [second_rolled_v]
                    second_pitched_v = second_rolled_v

            first_pitched_v = first_rolled_v

    # ===========================================

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
    branching_angle_roll = 0 #valid range [0,360]
    branching_angle_pitch = 20 #valid range [10,90]

    if (cost>threshold):
        #backtrack
        pass
    else:  #
        #todo - decide on the values of parameters by sampling the PDF
        output_point_list = generate_lsystem_tree(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch)
        cost = estimate_error(output_point_list, growth_space)
    return growth_param_list

def main():
    age = 4
    no_1st_ord_branches = 2
    no_2nd_ord_branches = 3
    branching_angle_roll = 30
    branching_angle_pitch = 21
    print generate_lsystem_tree(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch)

if __name__ == "__main__":
    main()
