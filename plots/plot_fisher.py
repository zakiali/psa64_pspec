#! /usr/bin/env python
import numpy as n, pylab as p, capo

npz = n.load('data/oqe.npz')
F = npz['FC']
kpl = npz['kpl']

fig = p.figure(figsize=(5,5))
#im1 = capo.arp.waterfall(F, drng=4, extent=(kpl[0],kpl[-1],kpl[-1],kpl[0]))
F = n.absolute(F)
F = n.ma.masked_less_equal(F,0)
F = n.ma.log10(F)
im1 = p.imshow(F, vmax=F.max(), vmin=F.max()-4, extent=(kpl[0],kpl[-1],kpl[-1],kpl[0]), aspect='auto', interpolation='nearest')
fig.subplots_adjust(left=.15, top=.95, bottom=.11, wspace=.2, hspace=.1, right=0.9)
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
p.ylabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
cbar_axis=fig.add_axes([.91, .15, .03, .75])
fig.colorbar(im1, cax=cbar_axis)

p.show()
