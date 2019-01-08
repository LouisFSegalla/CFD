import numpy as np
import pylab as plt
from glob import glob
import os
plt.ioff()


filenames = glob("./out/Solução no tempo*.txt")
filenames.sort()

for filename in filenames:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    f = open(filename, 'r')

    xCells = int(f.readline().split(":")[1])
    yCells = int(f.readline().split(":")[1])
    numGhostCells = int(f.readline().split(":")[1])
    time = float(f.readline().split(":")[1])
    f.close()

    x,y,u = np.loadtxt(filename,skiprows=4,unpack='true')
    x = np.reshape(x,(xCells + 2*numGhostCells, yCells + 2*numGhostCells))
    y = np.reshape(y,(xCells + 2*numGhostCells, yCells + 2*numGhostCells))
    u = np.reshape(u,(xCells + 2*numGhostCells, yCells + 2*numGhostCells))
    ax.set_aspect('equal', 'datalim')
    ax.contourf(x,y,u)
    ax.set_title("Time = %5.3f"%time)
    fig.savefig(filename.replace(".txt",".png"))


#os.system("eog " + filename[0].replace(".txt",".png"))
