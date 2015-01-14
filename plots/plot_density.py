#! /usr/bin/env python
import numpy as n, pylab as p
import os, gzip

#colors = 'kbmcrg'

densityfile = 'data/density.npz.gz'

f = gzip.open(densityfile)
fout = open(densityfile[:-len('.gz')], 'w')
s = f.read()
fout.write(s)
fout.close()
npz = n.load(densityfile[:-len('.gz')])
print npz.files
density = npz['density']
density[3500] = density[3501] # mask out point where flagged data accumulates
rng = npz['rng']
import capo as C
p.figure(figsize=(12,3))
C.arp.waterfall(density, extent=rng, origin='lower', cmap='gist_earth_r', mode='lin', mx=40)
p.colorbar(ticks=[0,10,20,30,40])
p.ylim(-500,500)
#p.ylim(-4500,4500)
p.xlim(0,n.pi)
#p.xlim(0,2)
p.grid()
p.subplots_adjust(right=.99, bottom=.2)
p.xlabel('LST [radians]')
p.ylabel('Visibility (real) [Jy]')
p.show()

for i in [500,1000,1500,2000]:
    p.plot(density[:,i] - n.average(density[n.where(density[:,i] > 0),i]))
p.show()
