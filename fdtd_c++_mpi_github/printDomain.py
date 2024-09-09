import csv
import numpy as np
import matplotlib.pyplot as plt


eps = []
with open("eps.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=",")
    for i,line in enumerate(reader):
        if i == 0:
            Lx, Ly, Nx, Ny = [float(s) for s in line]
        elif len(line) == 0:
            break
        else:
            eps.append([float(e) for e in line[:-1]])

mu = []
with open("mu.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=",")
    for i,line in enumerate(reader):
        if i == 0:
            Lx, Ly, Nx, Ny = [float(s) for s in line]
        elif len(line) == 0:
            break
        else:
            mu.append([float(u) for u in line[:-1]])

sigmaE = []
with open("sigmaE.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=",")
    for i,line in enumerate(reader):
        if i == 0:
            Lx, Ly, Nx, Ny = [float(s) for s in line]
        elif len(line) == 0:
            break
        else:
            sigmaE.append([float(u) for u in line[:-1]])

sigmaH = []
with open("sigmaH.csv") as csvfile:
    reader = csv.reader(csvfile, delimiter=",")
    for i,line in enumerate(reader):
        if i == 0:
            Lx, Ly, Nx, Ny = [float(s) for s in line]
        elif len(line) == 0:
            break
        else:
            sigmaH.append([float(u) for u in line[:-1]])


eps = np.array(eps)
mu = np.array(mu)
sigmaE = np.array(sigmaE)
sigmaH = np.array(sigmaH)


        
extent = [-1000*Lx/2, 1000*Lx/2, -1000*Ly/2, 1000*Ly/2]
fig, axs = plt.subplots(2,2, figsize=(12,8), sharex=True)

im1 = axs[0,0].imshow(np.rot90(eps,k=1), extent=extent)
axs[0,0].set_title("Epsilon", fontsize=18)
#axs[0,0].set_xlabel("x [nm]", fontsize=16)
axs[0,0].set_ylabel("y [nm]", fontsize=16)
axs[0,0].tick_params(axis='both', which='major', labelsize=14)
fig.colorbar(im1, ax=axs[0,0])

im2 = axs[1,0].imshow(np.rot90(mu,k=1), extent=extent)
axs[1,0].set_title("Mu", fontsize=18)
axs[1,0].set_xlabel("x [nm]", fontsize=16)
axs[1,0].set_ylabel("y [nm]", fontsize=16)
axs[1,0].tick_params(axis='both', which='major', labelsize=14)
fig.colorbar(im2, ax=axs[1,0])

im3 = axs[0,1].imshow(np.rot90(sigmaE,k=1), extent=extent)
axs[0,1].set_title("sigmaE", fontsize=18)
#axs[0,1].set_xlabel("x [nm]", fontsize=16)
axs[0,1].set_ylabel("y [nm]", fontsize=16)
axs[0,1].tick_params(axis='both', which='major', labelsize=14)
fig.colorbar(im3, ax=axs[0,1])

im4 = axs[1,1].imshow(np.rot90(sigmaH,k=1), extent=extent)
axs[1,1].set_title("sigmaH", fontsize=18)
axs[1,1].set_xlabel("x [nm]", fontsize=16)
axs[1,1].set_ylabel("y [nm]", fontsize=16)
axs[1,1].tick_params(axis='both', which='major', labelsize=14)
fig.colorbar(im4, ax=axs[1,1])

plt.savefig("Domain.png")
