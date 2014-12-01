#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C

d, kwds = a.img.from_fits('im0000.dbm.fits')
d = d[150:450,150:450]
#d = d[200:400,200:400]
d = a.img.recenter(d,(d.shape[0]/2,d.shape[1]/2))
#C.arp.waterfall(d, mode='log'); p.show()
_d = n.fft.fft2(d) 
# get rid of remaining n/s baselines
_d[:,:15] = 0
_d[:,-15:] = 0

#_d = n.where(n.abs(_d) > 1e-2, 1., 0)
ans = n.sum(n.abs(_d)**2) / n.sum(n.abs(_d))**2
print 1./ans

C.arp.waterfall(_d, mode='log')
p.show()

