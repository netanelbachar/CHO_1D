import numpy as np

from MD_QHO_Functions import *

coord = []
ar = np.linspace(1, 8, 8)
for i in range(len(ar)):
    coord.append(["P", str(ar[i]), str(float(0)), str(float(0))])

filename = "ringcarmi.xyz"
xyz = open(filename, 'w')
xyz.write(str(beads))
xyz.write("\n")
xyz.write("time=0.1")
xyz.write("\n")
for j in range(len(coord)):
    lst = coord[j]
    strlist = ' '.join(lst)
    xyz.write(strlist)
    xyz.write("\n")














