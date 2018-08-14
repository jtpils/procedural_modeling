# Created by Like on 7 Mar 2018

import numpy
from pyquaternion import Quaternion
import os
from converters.lpy2mtg2pts import lsystem_run

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

    return ngenerate_lsystem_tree_points(   age=int(params[0]),
                                            no_1st_ord_branches=int(params[1]),
                                            no_2nd_ord_branches=int(params[2]),
                                            branching_angle_roll=params[3],
                                            branching_angle_pitch=params[4] )

def ngenerate_lsystem_tree_points(  age=0,
                                    no_1st_ord_branches=2,
                                    no_2nd_ord_branches=2,
                                    branching_angle_roll=180.0,
                                    branching_angle_pitch=30.0 ):

    'Parse L-System string output into a list of points of cylinder base centers and their radii forming the trunk-branch skeleton representation of a tree'

    output_point_list = lsystem_run(age, no_1st_ord_branches, no_2nd_ord_branches, branching_angle_roll, branching_angle_pitch)
    # return parse_lstring(lstring)
    #parse l-system output string into 3d points matching Xfrog object size

    '''
    # temporary hardcoding of tree growth: age 1-2 trunk grows taller by 2 m/yr, age 3 first branch out by 1 m/yr, age 4 second branch out by 0.5 m/yr
    # ========================
    trunk_growth_rate = 2.0
    branch_growth_rate = 1.0
    pitch_angle = branching_angle_pitch/180.0 * 3.14159 #radian
    roll_angle = branching_angle_roll/180.0 * 3.14159 #radian

    current_post = numpy.array([0.0, 0.0, 0.0]) #tracking turtle's position
    current_dir = numpy.array([0.0, 0.0, 1.0])  #tracking turtle's orientation/heading (vector H)
    current_up = numpy.array([0.0, 1.0, 0.0])   #tracking turtle's up (vector U)
    current_left = numpy.array([1.0, 0.0, 0.0]) #tracking turtle's left (vector L)

    branching_stack = []

    if age == 0:
        return [current_post]

    output_point_list = [current_post]  #ground
    for t in range (1, age+1): # growing trunk
        current_post = current_post + current_dir*trunk_growth_rate   #do not use operator +=, as it overwrites all values into one current_post memory location
        output_point_list += [current_post]

    if age > 4: #1st order branches start growing
        #push state into stack
        branching_stack.append([current_post, current_dir, current_up, current_left])
        #print branching_stack

        pitch1 = Quaternion(axis=current_left,angle=pitch_angle)
        roll1 = Quaternion(axis=current_dir,angle=roll_angle)
        current_dir = pitch1.rotate(current_dir)
        current_up = pitch1.rotate(current_up)
        current_left = pitch1.rotate(current_left)
        first_branching_node = current_post

        for s in range (0, no_1st_ord_branches): #1st order branches
            current_dir = roll1.rotate(current_dir)
            current_up = roll1.rotate(current_up)
            current_left = roll1.rotate(current_left)

            for t in range(age-4): #growing 1st order branch
                current_post = first_branching_node + (t+1)*current_dir*branch_growth_rate
                output_point_list += [current_post]

            if age > 5: #2nd order branches start growing
                #push state into stack
                #print("2nd order - Current_post=%s" % (current_post))
                branching_stack.append([current_post, current_dir, current_up, current_left])
                #print branching_stack

                pitch2 = Quaternion(axis=current_left,angle=pitch_angle)
                roll2 = Quaternion(axis=current_dir,angle=roll_angle)
                current_dir = pitch2.rotate(current_dir)
                current_up = pitch2.rotate(current_up)
                current_left = pitch2.rotate(current_left)
                second_branching_node = current_post

                for u in range (0, no_2nd_ord_branches): #2nd order branch
                    current_dir = roll2.rotate(current_dir)
                    current_up = roll2.rotate(current_up)
                    current_left = roll2.rotate(current_left)

                    for v in range(age-5): #growing 2nd order branch
                        current_post = second_branching_node + (v+1)*current_dir*branch_growth_rate
                        output_point_list += [current_post]

                #pop stack to previous state
                [current_post, current_dir, current_up, current_left] = branching_stack.pop()
                #print branching_stack

        #pop stack to previous state
        [current_post, current_dir, current_up, current_left] = branching_stack.pop()
        #print branching_stack

    # ===========================================
    '''

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
    #output = generate_lsystem_tree_points(numpy.array([20,5,8,75.0,45.0]))
    output = ngenerate_lsystem_tree_points( age=2,
                                            no_1st_ord_branches=2,
                                            no_2nd_ord_branches=3,
                                            branching_angle_roll=75.0,
                                            branching_angle_pitch=45.0)
    print output
    file = open('output.obj', 'w')
    for item in output:
        file.write("v %d %d %d\n" % (item[0], item[1], item[2]))
    file.close()

if __name__ == "__main__":
    main()
