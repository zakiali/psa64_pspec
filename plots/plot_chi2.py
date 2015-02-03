#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import aipy as ap
import capo

lsts = np.loadtxt('data/omnical/psa6266_lsts.txt')
f = np.load('data/omnical/chisqdof_6266.npz')
fig, ax = plt.subplots()
#ax.imshow(f['arr_0'], vmin = 0, vmax = 3, extent=[100, 200, lsts[0], lsts[-1]], interpolation = 'none')
#ax.imshow(np.transpose(f['arr_0']), vmin = 0, vmax = 3, extent=[lsts[0], lsts[1],100,200], interpolation = 'none')
capo.arp.waterfall(np.transpose(f['arr_0']), mx=1, drng=1, extent=[lsts[0], lsts[1], 100, 200])
#ax.set_aspect(1)
plt.show()
