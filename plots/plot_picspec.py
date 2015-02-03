#!/usr/bin/env python
import numpy as n
import pylab as p

d = n.load('data/gianni/PicA_normalized_spectrum.npz')
freqs = d['freq']
flux = d['flux']
flux_err = d['flux_err']
pic = d['model']

p.semilogy(freqs,pic,'k', linewidth=3)
p.errorbar(freqs,flux,flux_err)
#p.plot(freqs,flux)
p.gca().set_yscale('log', nonposy='clip')
p.ylim(300,600)
p.xlabel('Frequency [MHz]', size='large')
p.ylabel('Source Fluc [Jy]', size='large')
p.yticks(n.arange(3e2,7e2,1e2), ['300','400','500','600'])
p.xlim(120,173)
p.show()
