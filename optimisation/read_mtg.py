import numpy as np

mtg_file = 'lstring_output.mtg'

print "Reading file: %s" % mtg_file

with open(mtg_file) as f:
    lines = f.readlines()

i=0
while(1):
    if lines[i].startswith("ENTITY-CODE"):
        break
    else:
        i=i+1
# Columns are 1=y, 2=radius, 3=z, 4=z
pts = np.genfromtxt(mtg_file,skip_header=i+1,usecols=(4,1,3))

print pts[0:5,:]
