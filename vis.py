#!/usr/bin/python3
import numpy as n
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation,writers
import glob

datfiles = glob.glob('output.*.dat')

data = dict()

for datfile in datfiles:
    t = datfile.split('.')[1]
    data[t] = n.loadtxt(datfile)

sortedkeys = sorted(data.keys())
initial = data[sortedkeys[0]]

fig = plt.figure(figsize=(15,10))
ax = plt.axes()



def animate(i):
    dat = data[sortedkeys[i]]
    ax.clear()
    ax.plot(dat[0], dat[1], '.')
    print (sortedkeys[i])

    rms = str(n.sqrt(n.mean(n.square(dat[1]))))

    plt.title(sortedkeys[i] + ' rms=' + rms)
    plt.grid()

    i += 1
    if i == len(sortedkeys):
        i = 0
Writer = writers['ffmpeg']
writer = Writer(fps = 10 )
ani = FuncAnimation(fig, animate,frames=range(len(sortedkeys)), interval=500)
ani.save('animation.mp4',writer=writer)
