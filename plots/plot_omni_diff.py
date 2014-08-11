import numpy as n, aipy as a
import pylab as p
import sys,capo,glob,os

#args = sys.argv[1:]
args = glob.glob('data/raw/omnical_outputs/zen.2456268.3*.uvcRREcACO')
aa = a.cal.get_aa('psa6240_v003', .1,.1,.1)
seps=['0,1']
sep2ij = {}
ij2sep = {}
s = ''
for sep in seps:
    sep2ij[sep] = capo.dfm.grid2ij(aa.ant_layout)[0][sep].split(',')
    s += capo.dfm.grid2ij(aa.ant_layout)[0][sep] 
    
    for ij in sep2ij[sep]:
        ij2sep[ij] = sep

conj = {}
toconj = capo.dfm.grid2ij(aa.ant_layout)[1]
for k in toconj.keys():
    conj['%d_%d'%a.miriad.bl2ij(k)] = toconj[k]

pol='xx'

times,data,flags = capo.arp.get_dict_of_uv_data(args, s, pol)
_d = n.zeros((len(times),203), dtype=n.complex64)
_w = n.zeros((len(times),203))

for bl in data.keys():
    i,j = a.miriad.bl2ij(bl)
    if conj['%d_%d'%(i,j)]:
        _d += n.conj(data[bl][pol])
    else:
        _d += data[bl][pol]
    _w += n.logical_not(flags[bl][pol])
_d = n.where(n.isnan(_d/_w), 0, _d/_w)


#get mdl data
chunk = args[0].split('.')[1] + '.' + args[0].split('.')[2][0] 
print glob.glob(os.path.dirname(args[0]) + '/data*%s*xx*.npz'%chunk)
npzfile = glob.glob(os.path.dirname(args[0]) + '/data*%s*xx*.npz'%chunk)
npz = n.load(npzfile[0])
sepdict = {'0,1':'sep24', '-1,1':'sep6', '1,1':'sep35'}

mdl_data = npz[sepdict[sep]]
mdl_data = n.where(n.isnan(mdl_data), 0, mdl_data)

#plot

p.figure(1)
p.subplot(131)
capo.arp.waterfall(mdl_data, mx=3, drng=3); p.colorbar(shrink=.5)
p.subplot(132)
capo.arp.waterfall(_d, mx=3, drng=3); p.colorbar(shrink=.5)
p.subplot(133)
capo.arp.waterfall(mdl_data - _d, mx=2); p.colorbar(shrink=.5)


p.figure(2)
p.subplot(131)
capo.arp.waterfall(mdl_data, mode='phs'); p.colorbar(shrink=.5)
p.subplot(132)
capo.arp.waterfall(_d, mode='phs'); p.colorbar(shrink=.5)
p.subplot(133)
capo.arp.waterfall((mdl_data - _d), mode='phs'); p.colorbar(shrink=.5)


p.figure(3)
p.subplot(131)
capo.arp.waterfall(mdl_data, mx=2000, drng=4000, mode='real'); p.colorbar(shrink=.5)
p.subplot(132)
capo.arp.waterfall(_d, mode='real', mx=2000, drng=4000); p.colorbar(shrink=.5)
p.subplot(133)
capo.arp.waterfall(mdl_data - _d, mode='real', mx=2000, drng=4000); p.colorbar(shrink=.5)


p.figure(4)
p.subplot(131)
capo.arp.waterfall(mdl_data, mode='imag', mx=2000, drng=4000); p.colorbar(shrink=.5)
p.subplot(132)
capo.arp.waterfall(_d, mode='imag', mx=2000, drng=4000); p.colorbar(shrink=.5)
p.subplot(133)
capo.arp.waterfall(mdl_data - _d, mode='imag', mx=2000, drng=4000); p.colorbar(shrink=.5)

#p.figure(5)
#p.subplot(131)
#capo.arp.waterfall(n.angle(mdl_data), mode='lin'); p.colorbar(shrink=.5)
#p.subplot(132)
#capo.arp.waterfall(n.angle(_d), mode='lin'); p.colorbar(shrink=.5)
#p.subplot(133)
#capo.arp.waterfall((n.angle(mdl_data) - n.angle(_d)), mode='phs'); p.colorbar(shrink=.5)



p.show()
