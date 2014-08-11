import capo
import numpy as n
import pylab as p

#delay filtered
d_xx = n.load('data/raw/omnical_lstbin_v2Bf_xx.npz')
d_yy = n.load('data/raw/omnical_lstbin_v2Bf_yy.npz')

sep = 'sep24' # this is the 0,1 spacing

freqs = n.arange(.1,.2,.1/203)
lsts = d_xx['lsts']*12/n.pi
extent=(freqs[0], freqs[-1], lsts[-1], lsts[0])
p.subplot(221)
p.title('XX')
capo.arp.waterfall(d_xx['sep24'], extent=extent, mx=2, drng=4)
p.subplot(222)
p.title('YY')
capo.arp.waterfall(d_yy['sep24'], extent=extent, mx=2, drng=4); p.colorbar(shrink=.5)
p.subplot(223)
capo.arp.waterfall(d_xx['sep24'], mode='phs', extent=extent)
p.subplot(224)
capo.arp.waterfall(d_yy['sep24'], mode='phs', extent=extent); p.colorbar(shrink=.5)

p.savefig('wbdf_waterfalls.png', format='png')
p.show()
