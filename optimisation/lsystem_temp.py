from pyquaternion import Quaternion
import numpy

def generate_lsystem_tree_points(p):
	'Parse L-System string output into a list of points of cylinder base centers and their radii forming the trunk-branch skeleton representation of a tree'
	age = int(p[0]);
	no_1st_ord_branches = int(p[1]);
	no_2nd_ord_branches = int(p[2]);
	branching_angle_roll = p[3];
	branching_angle_pitch = p[4];
	trunk_growth_rate = p[5]
	branch_growth_rate = p[6]
    # temporary hardcoding of tree growth: age 1-2 trunk grows taller by 2 m/yr, age 3 first branch out by 1 m/yr, age 4 second branch out by 0.5 m/yr
    # ========================
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

    # ===========================================

        return output_point_list
