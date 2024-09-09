import csv
import numpy as np
import matplotlib.pyplot as plt
import os

CSV = "I_planeWave.csv"
folder = "I_planeWave"

if not os.path.exists(folder):
             os.makedirs(folder)

#I_max = 0
#I_line = []
#with open(CSV) as csvfile:
#    reader = csv.reader(csvfile, delimiter=",")
#    for i,line in enumerate(reader):
#        if len(line) != 2:
#            I_line = [float(l) for l in line[:-1]]
#            for I in I_line:
#                if I>I_max:
#                    I_max = I


I_snap = []
with open(CSV) as csvfile:
    reader = csv.reader(csvfile, delimiter=",")
    for i,line in enumerate(reader):
        if i == 0:
            Lx, Ly = [float(l) for l in line]
        elif len(line) == 2:
            if len(I_snap) != 0:
                I = np.array(I_snap)#/I_max
                plt.imshow(np.rot90(I), extent=[-Lx/2, Lx/2, -Ly/2, Ly/2])
                plt.title(t)
                #plt.clim(0,1)
                plt.colorbar()
                plt.savefig(folder+"/"+str(index)+".png")
                plt.clf()
            index, t = [float(l) for l in line]
            I_snap = []
        else:
            I_snap.append([float(l) for l in line[:-1]])