#!/usr/bin/env python
import numpy as n
#import matplotlib.pyplot as plt
import pylab as p
import aipy as a
import capo

lsts = n.loadtxt('data/omnical/psa6266_lsts.txt')
startlst = (lsts[0] - 2*n.pi) * 12/n.pi
endlst = lsts[-1]*12/n.pi

f = n.load('data/omnical/chisqdof_6266.npz')
fig = p.figure(figsize=(7,4))
ax1 = p.subplot(1,1,1)
#ax.imshow(f['arr_0'], vmin = 0, vmax = 3, extent=[100, 200, lsts[0], lsts[-1]], interpolation = 'none')
#ax.imshow(n.transpose(f['arr_0']), vmin = 0, vmax = 3, extent=[lsts[0], lsts[1],100,200], interpolation = 'none')
print n.max(n.where(f['arr_0'])), n.min(n.where(f['arr_0']))
print n.mean(n.where(f['arr_0']))
#capo.arp.waterfall(n.transpose(f['arr_0']), mx=300, drng=100,  extent=[lsts[0], lsts[1], 100, 200], mode='lin'); p.colorbar(shrink=.5)
#im1 = capo.arp.waterfall(n.transpose(f['arr_0']), mx=.5, drng=.5, extent=[startlst,endlst, 100, 200], mode='log'); p.colorbar(shrink=.5)
im1 = p.imshow(n.log10(n.abs(n.transpose(f['arr_0']))), vmax=.5, vmin=0, extent=[startlst,endlst, 100, 200], interpolation='nearest', aspect='auto')
p.xlabel('LST [ Hours ]')
p.ylabel('Frequency [ MHz ]')

fig.subplots_adjust(left=.1, top=.95, bottom=.15, wspace=.2, hspace=.1, right=0.8)
cbar = fig.add_axes([.85,.2,.045,.7])
fig.colorbar(im1, cax=cbar)
#cbar_ax1 = fig.add_axes([.9,.1,.05, .3])
#fig.colorbar(im1, cax=cbar_ax1)
#ax.imshow(n.transpose(f['arr_0']), vmin = 0, vmax = 3, extent=[lsts[0], lsts[1],100,200], interpolation = 'none')
#ax.set_aspect(1)
p.show()
