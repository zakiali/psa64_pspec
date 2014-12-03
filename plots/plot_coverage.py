#! /usr/bin/env python
"""
Script for displaying a projection of a spherical (Healpix) data set stored
in a *.fits file.

Author: Aaron Parsons
"""

import aipy as a, numpy as n, pylab as p, sys, os, ephem, optparse

# Try to import basemap module, but on failure use the above class
from mpl_toolkits.basemap import Basemap



o = optparse.OptionParser()
o.set_usage('plot_map.py [options] mapfile')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, src=True, 
    cmap=True, max=True, drng=True)
o.add_option('-p', '--projection', dest='projection', default='moll',
    help='Map projection to use: moll (default), mill, cyl, robin, sinu.')
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plotting mode, can be log (default), lin.')
o.add_option('--interpolation', dest='interpolation', default='nearest',
    help='Sub-pixel interpolation.  Can be "nearest" or "bicubic".  Default nearest.')
o.add_option('--iepoch', dest='iepoch', type='float', default=ephem.J2000,
    help='Epoch of input coordinates (in map).  Default J2000.')
o.add_option('--oepoch', dest='oepoch', type='float', default=ephem.J2000,
    help='Epoch of output coordinates (plotted).  Default J2000.')
o.add_option('--res', dest='res', type='float', default=0.25,
    help="Resolution of plot (in degrees).  Default 0.25.")
o.add_option('--nside', dest='nside', type='int',
    help="Manually set NSIDE (possibly degrading map) to a power of 2.")
opts,args = o.parse_args(sys.argv[1:])
args = ['data/dOC-150MHz.fits', 'data/bm_coverage.fits']

if not os.path.exists(args[1]):
    npz = n.load('data/coverage.npz')
    lsts = npz['times']
    coverage = n.max(npz['cnt'], axis=1)
    aa = a.cal.get_aa('psa6240_v003', .001, .15, 1)
    nh = a.healpix.HealpixMap(nside=32)
    eq = nh.px2crd(n.arange(nh.npix()), ncrd=3)
    for lst,cnt in zip(lsts,coverage):
        print lst
        m = a.coord.eq2top_m(-lst, aa.lat)
        tx,ty,tz = n.dot(m, eq)
        bmx = aa[0].bm_response((tx,ty,tz), pol='x')**2
        bmy = aa[0].bm_response((tx,ty,tz), pol='y')**2
        bm = n.where(tz > 0, 0.5 * (bmx + bmy), 0)
        #nh.map += cnt * bm.flatten()
        nh.map += bm.flatten()
    nh.to_fits(args[1])
print 'Reading %s' % args[1]
nh = a.healpix.HealpixMap(fromfits=args[1])

p.figure(figsize=(10,5))
cmap = p.get_cmap(opts.cmap)
map = Basemap(projection=opts.projection,lat_0=0,lon_0=270, rsphere=1.)
lons,lats,x,y = map.makegrid(360/opts.res,180/opts.res, returnxy=True)
# Mask off parts of the image to be plotted that are outside of the map
lt = lats[:,0]
ln1 = n.ones_like(lt) * (lons[lons.shape[0]/2,0])
ln2 = n.ones_like(lt) * (lons[lons.shape[0]/2,-1])
x1,y1 = map(ln1,lt); x2,y2 = map(ln2,lt)
x = n.ma.array(x)
for c,(i,j) in enumerate(zip(x1,x2)): x[c] = n.ma.masked_outside(x[c], i, j)
mask = x.mask
lons = 360 - lons
lats *= a.img.deg2rad; lons *= a.img.deg2rad
print 'Reading %s' % args[0]
h = a.map.Map(fromfits=args[0])
print 'SCHEME:', h.scheme()
print 'NSIDE:', h.nside()
h.set_interpol(opts.interpolation != 'nearest')

crd = a.coord.radec2eq(n.array([lons.flatten(), lats.flatten()]))
m = a.coord.convert_m('eq', 'ga', iepoch=opts.oepoch, oepoch=opts.iepoch)
x,y,z = n.dot(m, crd)
gsm = h[x,y,z]
gsm.shape = lats.shape
m = a.coord.convert_m('eq', 'eq', iepoch=opts.oepoch, oepoch=opts.iepoch)
x,y,z = n.dot(m, crd)
nh.set_interpol(True)
coverage = nh[x,y,z]
coverage.shape = lats.shape


# Generate map grid/outline
map.drawmapboundary()
map.drawmeridians(n.arange(-180, 180, 30))
map.drawparallels(n.arange(-90,90,30)[1:], labels=[0,1,0,0], labelstyle='+/-')
# Set up data to plot
gsm = n.log10(n.abs(gsm))
max,min = 4,2
gsm = gsm.clip(min, max)
gsm = n.ma.array(gsm, mask=mask)
coverage /= n.max(coverage)
#coverage = n.log10(coverage/n.max(coverage))
coverage = n.ma.array(coverage, mask=mask)
print coverage.max()
map.imshow(gsm, vmax=max, vmin=min, cmap=cmap, interpolation=opts.interpolation)
#map.imshow(coverage, cmap='gray', interpolation=opts.interpolation, alpha=0.5)
lons,lats,x,y = map.makegrid(360/opts.res,180/opts.res, returnxy=True)
#map.contour(x,y, coverage, n.arange(.2,1,.2), cmap='hot', interpolation=opts.interpolation)
manual = [(3.02,1.08), (3.37,1.20), (2.15,1.13), (2.20,1.50), (2.11,2.09)]
cs = map.contour(x,y, coverage, [1e-2, 1e-1, .25, .5, .75], colors='w')
p.clabel(cs, inline=1, fontsize=10, fmt='%4.2f', manual=manual)

p.subplots_adjust(.05,.05,.95,.95)

p.show()
