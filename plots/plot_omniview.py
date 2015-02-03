#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import aipy as ap
f = np.load('data/omnical/omniview.npz')
fig, axes = plt.subplots(nrows=1, ncols=len(f['arr_0']), sharey=True, sharex=True, figsize=(10,6))
for plot_i, datagrps in enumerate(f['arr_0']):

    #plt.subplot('12%i'%(2-plot_i))
    for (data, marker, color, blvec) in datagrps:
        axes[plot_i].scatter(np.real(data), np.imag(data), marker = marker, color=color)
        axes[plot_i].set(adjustable='datalim', aspect=1)#set_aspect(1)
        axes[plot_i].grid(True)
plt.axis(8000*np.array([-1, 1, -1, 1]))
#plt.xlabel('Real [Jy]', size='Large')
#plt.ylabel('Imaginary [Jy]', size='Large')
axes[0].set_ylabel('Imaginary [Jy]', size='large')
axes[0].set_xlabel('Real [Jy]', size='large')
axes[0].set_title('Rough Calibration',size='large')
#axes[1].set_ylabel('Imaginary [Jy]', size='large')
axes[1].set_xlabel('Real [Jy]', size='large')
axes[1].set_title('Omnical',size='large')

plt.show()

