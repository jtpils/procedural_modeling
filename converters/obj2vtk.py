#29 March 2018
#Convert .obj to .vtk for all relevant files in one folder, for data points only


import os

input_path = "C:/Users/gobeawanl/Documents/VS1 Tree Modelling/Work/Daniel/20180328 Optimisation 5 params 2000 init conds 122min runtime/objs/"
output_path = "C:/Users/gobeawanl/Documents/VS1 Tree Modelling/Work/Daniel/20180328 Optimisation 5 params 2000 init conds 122min runtime/vtks/"
file_names = [file_name for file_name in os.listdir(input_path) if os.path.isfile(os.path.join(input_path, file_name))]
vtktype = "POLYDATA"    #type of data in VTK
datatype = "integer"    #datatype of numbers in VTK

def find_no_points(file_name):
    with open(file_name) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

for file_name in file_names:
    fi = open(os.path.join(input_path, file_name), "r")
    file_out = file_name.rsplit('.', 1)
    fo = open(os.path.join(output_path, file_out[0]+".vtk"), "w")

    no_points = find_no_points(os.path.join(input_path, file_name))
    fo.write("# vtk DataFile Version 3.0\nTree points\nASCII\nDATASET "+vtktype+"\n\nPOINTS "+str(no_points)+" "+datatype+" ")
    for aPoint in fi:
        fo.write("%s" % aPoint[2:])

    fi.close()
    fo.close()
