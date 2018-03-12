# Created by Like on 7 Mar 2018
# Last modified 12 Mar 2018

from pyquaternion import Quaternion

def estimate_age():
    age = 4
    return age

def lsystem_run(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch):
    'Pass known parameter values into L-system rules to produce Lstring output'
    pass #todo - Peng
    return lstring_output

def generate_lsystem_tree(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch):
    'Translate Lstring output into a list of points of cylinder base centers and their radii forming the trunk-branch representation of a tree'

    if age == 0:
        return [[0,0,0]]

    out = lsystem_run(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch)

    #parse l-system output string
    #todo
    #to match with Xfrog object size

    output_point_list = [[0,0,0]]
    for t in range (1, age+1): #trunk
        output_point_list = output_point_list + [[0,0,t]]

    v = [0,0,1]
    pitch_axis = [-1,0,0]
    roll_axis = [0,0,1]
    pitch_angle = branching_angle_pitch/180 * 3.14 #radian
    roll_angle = branching_angle_roll/180 * 3.14 #radian
    pitched_v = Quaternion(axis=pitch_axis,angle=pitch_angle).rotate(v)
    if age > 2
        for t in range (3, age+1): #1st order branch
            rolled_v = Quaternion(axis=roll_axis,angle=roll_angle).rotate(pitched_v)
            pitched_v = rolled_v
            output_point_list = output_point_list + [rolled_v]

    if age > 3
        for t in range (4, age+1): #2nd order branch
            output_point_list = output_point_list + []




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

    if (cost>threshold)
        #backtrack
        pass
    else  #
        #todo - decide on the values of parameters by sampling the PDF
        output_point_list = generate_lsystem_tree(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch)
        cost = estimate_error(output_point_list, growth_space)
    return growth_param_list
