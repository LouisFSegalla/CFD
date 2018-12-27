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
    numGhostCells = int(f.readline().split(":")[1])
    time = float(f.readline().split(":")[1])
    f.close()

    x,u = np.loadtxt(filename,skiprows=3,unpack='true')
    p, = ax.plot(x,u,lw=3)
    ax.set_xlim(-1,1)
    ax.set_ylim(-1.5,1.5)
    ax.set_title("Time = %5.3f"%time)
    fig.savefig(filename.replace(".txt",".png"))


    #os.system("eog " + filename[0].replace(".txt",".png"))
