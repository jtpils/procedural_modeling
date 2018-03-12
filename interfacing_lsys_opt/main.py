# Created by Like on 7 Mar 2018
# Last modified 12 Mar 2018

def estimate_age():
    age = 4
    return age

def lsystem_run(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch):
    'Pass known parameter values into L-system rules to produce Lstring output'
    pass #todo - Peng
    return lstring_output

def generate_lsystem_tree(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch):
    'Translate Lstring output into a list of points of cylinder base centers and their radii forming the trunk-branch representation of a tree'

    out = lsystem_run(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch)

    #parse l-system output string
    #todo
    #to match with Xfrog object size


    return output_point_list

def estimate_error(output_point_list, growth_space):

    return cost;

def optimise():
    age = estimate_age()
    cost = 0

    if (cost>threshold)
        #todo
        pass
    else
        #todo - decide on the values of parameters by sampling the PDF
        output_point_list = generate_lsystem_tree(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch)
        cost = estimate_error(output_point_list, growth_space)
    return growth_param_list
