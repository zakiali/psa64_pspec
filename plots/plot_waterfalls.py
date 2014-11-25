#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import glob

def get_data(filenames, antstr, polstr, verbose=False):
    lsts, dat, flg = [], {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            lst = uv['lst']
            if len(lsts) == 0 or lst != lsts[-1]: lsts.append(lst)
            bl = a.miriad.ij2bl(i,j)
            if not dat.has_key(bl): dat[bl],flg[bl] = [],[]
            dat[bl].append(d)
            flg[bl].append(f)
    return n.array(lsts), dat, flg

uvA_fg = glob.glob('data/lstdata/uvA_fg/*uvAA')
uvAL_fg = glob.glob('data/lstdata/uvA_fg/*uvALA')
uvA = glob.glob('data/lstdata/uvA/*uvAA')
uvAL = glob.glob('data/lstdata/uvA/*uvALA')
uvALS = glob.glob('data/lstdata/uvA/*uvALSA')

bl = a.miriad.ij2bl(0,26)
fig = p.figure(figsize=(10,5))
ax = None
for cnt, dset in enumerate([uvA_fg, uvAL_fg, uvA, uvAL, uvALS]):
    if ax is None:
        ax = p.subplot(1,5,cnt+1)
        p.ylabel('LST [hours]')
    else:
        ax2 = p.subplot(1,5,cnt+1, sharey=ax)
        p.setp(ax2.get_yticklabels(), visible=False)
    lsts,d,f = get_data(dset, 'cross', 'I', verbose=True)
    lsts = n.where(lsts > 5, lsts-2*n.pi, lsts)
    lsts = 12 * lsts/n.pi
    im = C.arp.waterfall(d[bl], mode='log', mx=4, drng=5, extent=(100,200,lsts[-1],lsts[0]))
    p.plot([100,200],[0,0],'k--')
    p.plot([100,200],[9.5,9.5],'k--')
    p.ylim(17,-2)
    p.xlim(110,185)
    p.xticks([125, 150, 175])
    if cnt == 2: p.xlabel('Frequency [MHz]')

fig.subplots_adjust(left=.08, top=.95, bottom=.12, wspace=.3, right=0.85)
cbar_ax = fig.add_axes([0.90, 0.15, 0.03, 0.7])
fig.colorbar(im, cax=cbar_ax)
p.xlabel(r'${\rm log}_{10}({\rm Jy})$', fontsize=14)

p.show()
