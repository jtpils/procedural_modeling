#29 March 2018
#Convert .obj to .vtk for all relevant files in one folder, for data points only


import os

path_name = "C:/Users/gobeawanl/Documents/VS1 Tree Modelling/Work/Daniel/20180328 Optimisation 5 params 2000 init conds 122min runtime/vtks/"
file_names = [file_name in os.listdir(path_name) if os.path.isfile(os.path.join(path_name, file_name))]
vtktype = "POLYDATA"    #type of data in VTK
datatype = "integer"    #datatype of numbers in VTK

for file_name in file_names:
    file fi = open(file_name, "r")
    file_out = file_name.rsplit('.', 1)
    file fo = open(file_out[0]+".vtk", "w")

    no_lines =

    fo.write("# vtk DataFile Version 3.0\nTree points\nASCII\nDATASET "+vtktype+"\n\nPOINTS "+no_lines+" "+datatype+" ")
    for :
        fi.read()
        fo.write("%s")

    fi.close()
    fo.close()
