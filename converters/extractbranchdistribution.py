import os
import math

max_depth_to_process = 1    #max=1 means only depth 0 and 1 will be processed
current_path = os.path.dirname(os.path.abspath(__file__))
input_file = os.path.abspath(os.path.realpath(os.path.join(current_path, 'mtg/20190222114041_hopea.mtg')))
#output_file = os.path.join(current_path, 'mtg/test_branchdistribution.txt')

content_list = []
with open(input_file, 'rt') as in_file:
    for line in in_file:
        if line.startswith(('/', '+', '^', '\t')):
            content_list.append(line.rstrip('\n'))
    #print content_list

prev_depth = cur_depth = 0
#branch_length_stack = []
last_point_in_segment_stack = []
accum_distance = 0
prev_x = prev_y = prev_z = 0.0
for element in content_list:
    #print element
    pieces = element.split('\t')
    cur_x = float(pieces[len(pieces)-3])
    cur_y = float(pieces[len(pieces)-4])
    cur_z = float(pieces[len(pieces)-1])
    #print 'cur x,y,z=', cur_x, ',', cur_y, ',', cur_z
    if '+' in element:
        #print 'DEBUG +'
        cur_depth = len(element[:element.index('+')])
        #print 'cur_depth=', cur_depth,
        if cur_depth > max_depth_to_process:
            continue
        if cur_depth > prev_depth:  #branch out from current segment. If invalid syntax where cur_depth > prev_depth+1, treat it as cur_depth = prev_depth+1
            print accum_distance, '[',
            #branch_length_stack.append(accum_distance)
            last_point_in_segment_stack.append([prev_x, prev_y, prev_z])
            accum_distance = math.sqrt(math.pow(cur_x-prev_x, 2) + math.pow(cur_y-prev_y, 2) + math.pow(cur_z-prev_z, 2))
        elif cur_depth <= prev_depth:     #end of previous branch(es), branch out from previous segment
            print accum_distance,
            if cur_depth > 0:
                for i in range(prev_depth-cur_depth):
                    #print branch_length_stack.pop(), ']',
                    last_point_in_segment_stack.pop()
                    print ']',
                print '] [',
                accum_distance = math.sqrt(math.pow(cur_x-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][0], 2)
                                           + math.pow(cur_y-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][1], 2)
                                           + math.pow(cur_z-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][2], 2))
            else:   #invalid branching cur_depth, assume ^<I instead of +I
                cur_depth = 0
                for i in range(prev_depth-1):
                    last_point_in_segment_stack.pop()
                    print ']',
                accum_distance = math.sqrt(math.pow(cur_x-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][0], 2)
                                           + math.pow(cur_y-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][1], 2)
                                           + math.pow(cur_z-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][2], 2))
                last_point_in_segment_stack.pop()
                print ']',

    elif '^' in element:
        #print 'DEBUG ^'
        cur_depth = len(element[:element.index('^')])
        #print 'cur_depth=', cur_depth,
        if cur_depth > max_depth_to_process:
            continue
        if cur_depth == prev_depth:
            accum_distance += math.sqrt(math.pow(cur_x-prev_x, 2) + math.pow(cur_y-prev_y, 2) + math.pow(cur_z-prev_z, 2))
        elif cur_depth > prev_depth: #treat ^<I as +I
            print accum_distance, '[',
            #branch_length_stack.append(accum_distance)
            last_point_in_segment_stack.append([prev_x, prev_y, prev_z])
            accum_distance = math.sqrt(math.pow(cur_x-prev_x, 2) + math.pow(cur_y-prev_y, 2) + math.pow(cur_z-prev_z, 2))
        else:   #cur_depth < prev_depth
            #store previous accum_distance for finished branches
            #print 'DEBUG cur_depth=', cur_depth, 'prev_depth=', prev_depth,
            print accum_distance,
            for i in range(prev_depth-cur_depth-1):
                #print branch_length_stack.pop(), ']',
                last_point_in_segment_stack.pop()
                print ']',
            #print 'DEBUG: len(branch_length_stack)=', len(branch_length_stack)
            if len(last_point_in_segment_stack) < 1:
                print 'DEBUG len(last_point_in_segment_stack)=', len(last_point_in_segment_stack)
            accum_distance = math.sqrt(math.pow(cur_x-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][0], 2)
                                       + math.pow(cur_y-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][1], 2)
                                       + math.pow(cur_z-last_point_in_segment_stack[len(last_point_in_segment_stack)-1][2], 2))
            last_point_in_segment_stack.pop()
            print ']',

    elif '/' in element:
        #print 'DEBUG /'
        cur_depth = len(element[:element.index('/')])
        #print 'cur_depth=', cur_depth,
        if cur_depth > max_depth_to_process:
            continue
        accum_distance = 0
        print '[',

    prev_depth = cur_depth
    prev_x = cur_x
    prev_y = cur_y
    prev_z = cur_z

print accum_distance, ']'
