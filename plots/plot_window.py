#! /usr/bin/env python
import numpy as n, pylab as p, capo
from matplotlib import gridspec

npz = n.load('data/oqe.npz')
W = npz['WC']
kpl = npz['kpl']

gs = gridspec.GridSpec(2, 2, width_ratios=[1,.05], height_ratios=[3,1])
fig = p.figure(figsize=(6,6))
ax = p.subplot(gs[0])
im = capo.arp.waterfall(W, drng=4, extent=(kpl[0],kpl[-1],kpl[-1],kpl[0]))
p.plot([kpl[0],kpl[-1]], [kpl[6], kpl[6]], 'k--')
p.setp(ax.get_xticklabels(), visible=False)
p.ylabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
p.xlim(kpl[0],kpl[-1])
p.ylim(kpl[-1],kpl[0])
ax = p.subplot(gs[1])
fig.colorbar(im, cax=ax)

p.subplot(gs[2])
p.semilogy(kpl, n.abs(W[6]))
fig.subplots_adjust(left=.18, top=.95, bottom=.12, wspace=.1, hspace=.1, right=0.9)
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
p.xlim(kpl[0],kpl[-1])
p.yticks([1e0, 1e-4, 1e-8, 1e-12, 1e-16])
p.grid(True)
p.show()
