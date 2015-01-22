import numpy as n
import pylab as p



normal_noise = n.random.normal(size=100000)
median = n.median(normal_noise)
res = n.abs(normal_noise - median)
sigma = n.median(res)
print median
print sigma
filtered = n.where(res>3*sigma, 0, normal_noise)

varfilt = n.var(filtered)
ovar = n.var(normal_noise)
print ovar, varfilt 
print 'correction= ', ovar/varfilt

#p.plot(filtered)
#p.plot(normal_noise)
#p.show()
#
#ohist, obins = n.histogram(normal_noise, bins=100)
#hist, bins = n.histogram(filtered, bins=100)
#p.plot(obins[:-1], ohist)
#p.plot(bins[:-1], hist)
#p.vlines(median, 0, 50)
#p.show()

