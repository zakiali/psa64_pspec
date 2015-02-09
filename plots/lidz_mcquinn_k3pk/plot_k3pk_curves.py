#! /usr/bin/env python
import numpy as n, pylab as p
import capo as C
import sys, re
import scipy.interpolate

def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK

re_z = re.compile(r'power_21cm_z(\d+\.\d+)\.dat')
#(0:02; 11:46), (0.15, 8.76 ), (0.21, 8.34), (0.54, 7.32), (0.71, 7.03), (0.82, 6.90), and (0.96, 6.77)

limits = {}
#limits['.05-.1'] = (n.arange(.05, .1, .001), 6100.,'c')
limits['.1-.2'] = (n.arange(.1, .2, .001), 3700.,'m')
limits['.2-.4'] = (n.arange(.2, .4, .001), 5680.,'b')

dat = {}
for filename in sys.argv[1:]:
    print 'Reading', filename
    d = n.array([map(float, L.split()) for L in open(filename).readlines()])
    ks, pk = d[:,0], d[:,1]
    z_file = float(re_z.match(filename).groups()[0])
    z = C.pspec.f2z(.160)
    k3pk = ks**3 / (2*n.pi**2) * pk# * mean_temp(z)**2
    #dat[filename] = (ks, k3pk)
    #dat[filename] = (ks, pk)
    #dat['%5.3f'%z_file] = scipy.interpolate.interp1d(ks, pk, kind='linear')
    dat[z_file] = scipy.interpolate.interp1d(ks, pk, kind='linear')
    #p.loglog(ks, k3pk * mean_temp(z)**2, label='%s,%f'%(filename,z))
    p.loglog(ks, k3pk * mean_temp(z)**2, label='%5.3f->%5.3f'%(z_file,z))
p.legend(loc='best')
p.show()

#zs = dat.keys(); zs.sort()
colors = 'kbgrcmy'
#fqs = n.arange(.150, .190, .01)
fq = .160
#color = colors[i%len(colors)]
Tbs = {}
flist = dat.keys(); flist.sort()
for f in flist:
    print f
    for L in limits:
        ks, lim, clr = limits[L]
        pk = n.average(dat[f](ks))
        k = n.average(ks)
        k3pk = k**3 / (2*n.pi**2) * pk
        Tb = n.sqrt(lim / k3pk)
        Tbs[L] = Tbs.get(L,[]) + [Tb]
    #p.loglog(dat[f][0], dat[f][1] * mean_temp(z)**2, color, label='%s,%f'%(f,z))

for L in Tbs:
    print Tbs[L]
    p.plot(Tbs[L])
p.show()
