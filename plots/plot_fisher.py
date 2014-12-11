#! /usr/bin/env python
import numpy as n, pylab as p, capo

npz = n.load('data/oqe.npz')
F = npz['FC']
kpl = npz['kpl']

fig = p.figure(figsize=(5,5))
capo.arp.waterfall(F, drng=4, extent=(kpl[0],kpl[-1],kpl[-1],kpl[0]))
fig.subplots_adjust(left=.18, top=.95, bottom=.12, wspace=.2, hspace=.1, right=0.95)
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
p.ylabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
p.show()
