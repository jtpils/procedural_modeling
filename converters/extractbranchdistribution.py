import os
import math

current_path = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.abspath(os.path.realpath(os.path.join(current_path, 'mtg/test2.mtg')))
output_file = os.path.join(current_path, 'mtg/test_branchdistribution.txt')

content_list = []
with open(input_file, 'rt') as in_file:
    for line in in_file:
        if line[0] is '/' or line[0] is '+' or line[0] is '^' or line[0] is '\t':
            content_list.append(line.rstrip('\n'))

prev_depth = cur_depth = 0
branch_length_stack = []
last_point_in_segment_stack = []
accum_distance = 0
prev_x = prev_y = prev_z = 0.0
for element in content_list:
    pieces = element.split('\t')
    cur_x = float(pieces[len(pieces)-3])
    cur_y = float(pieces[len(pieces)-4])
    cur_z = float(pieces[len(pieces)-1])
    print 'cur x,y,z=', cur_x, ',', cur_y, ',', cur_z
    if '+' in element:
        cur_depth = len(element[:element.index('+')])
        #print 'cur_depth=', cur_depth,
        if cur_depth > prev_depth:  #branch out from current segment. If invalid syntax where cur_depth > prev_depth+1, treat it as cur_depth = prev_depth+1
            branch_length_stack.append(accum_distance)
            last_point_in_segment_stack.append([prev_x, prev_y, prev_z])
            accum_distance = math.sqrt(math.pow(cur_x-prev_x, 2) + math.pow(cur_y-prev_y, 2) + math.pow(cur_z-prev_z, 2))
        elif cur_depth <= prev_depth:     #end of previous branch(es), branch out from previous segment
            if cur_depth < 0:  #incorrect depth
                cur_depth = 0
            for i in range(prev_depth-cur_depth):
                print branch_length_stack.pop(), ']',
                last_point_in_segment_stack.pop()
            if cur_depth == prev_depth:
                print accum_distance, '['
            accum_distance = math.sqrt(math.pow(cur_x-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][0], 2)
                                       + math.pow(cur_y-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][1], 2)
                                       + math.pow(cur_z-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][2], 2))

    elif '^<' in element:
        cur_depth = len(element[:element.index('^')])
        #print 'cur_depth=', cur_depth,
        if cur_depth >= prev_depth:
            accum_distance += math.sqrt(math.pow(cur_x-prev_x, 2) + math.pow(cur_y-prev_y, 2) + math.pow(cur_z-prev_z, 2))
        else:   #cur_depth < prev_depth
            #store previous accum_distance for finished branches
            for i in range(prev_depth-cur_depth):
                print branch_length_stack.pop(), ']'
                last_point_in_segment_stack.pop()
            accum_distance = branch_length_stack[len(branch_length_stack)-1] + math.sqrt(math.pow(cur_x-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][0], 2)
                                                                                         + math.pow(cur_y-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][1], 2)
                                                                                         + math.pow(cur_z-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][2], 2))
    elif '/' in element:
        cur_depth = 0 #len(element[:element.index('/')])
        #print 'cur_depth=', cur_depth,
        accum_distance = 0
        print '[', accum_distance

    prev_depth = cur_depth
    prev_x = cur_x
    prev_y = cur_y
    prev_z = cur_z
