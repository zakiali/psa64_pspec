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
#uvA = glob.glob('data/lstdata/uvA/*uvAA')
#uvAL = glob.glob('data/lstdata/uvA/*uvALA')
#uvALS = glob.glob('data/lstdata/uvA/*uvALSA')
uvA = glob.glob('data/lstdata/uvA_v2/*uvA')
uvAL = glob.glob('data/lstdata/uvA_v2/*uvGLA')
uvALS = glob.glob('data/lstdata/uvA_v2/*uvGLSA')

bl = a.miriad.ij2bl(0,26)
fig = p.figure(figsize=(10,8))
titles = ['LST Binned', 'Fringe Rate Filt.', 'Foreground Filt.', 'Fringe Rate Filt.', 'Baseline Avg.']
for cnt, (dset,title) in enumerate(zip([uvA_fg, uvAL_fg, uvA, uvAL, uvALS],titles)):
    #if cnt == 0:
    #    ax2 = p.subplot(2,5,1+5)
    #    ax1 = p.subplot(2,5,1, sharex=ax2)
    #else:
    #    axa = p.subplot(2,5,cnt+1, sharey=ax1)
    #    p.setp(axa.get_xticklabels(), visible=False)
    #    axb = p.subplot(2,5,cnt+1+5, sharey=ax2)
    #    p.setp(axa.get_yticklabels(), visible=False)
    #    p.setp(axb.get_yticklabels(), visible=False)
    lsts,d,f = get_data(dset, 'cross', 'I', verbose=True)
    print d.keys()
    d = n.array(d[bl])
    lsts = n.where(lsts > 5, lsts-2*n.pi, lsts)
    lsts = 12 * lsts/n.pi

    ax1 = p.subplot(2,5,cnt+1)
    d[:,76:78] = 0 # make sure orbcomm is white
    im1 = C.arp.waterfall(d, mode='log', mx=4, drng=5, extent=(100,200,lsts[-1],lsts[0]))
    p.setp(ax1.get_xticklabels(), visible=False)
    p.title(title, size=13)
    p.plot([100,200],[0,0],'k--')
    p.plot([100,200],[9.5,9.5],'k--')
    p.ylim(17,-2)
    p.xlim(110,185)
    p.xticks([125, 150, 175])
    if cnt == 0: p.ylabel('LST [hours]')
    else: p.setp(ax1.get_yticklabels(), visible=False)

    ax2 = p.subplot(2,5,cnt+6)
    im2 = C.arp.waterfall(d, mode='phs', mx=n.pi, drng=2*n.pi, extent=(100,200,lsts[-1],lsts[0]))
    p.plot([100,200],[0,0],'k--')
    p.plot([100,200],[9.5,9.5],'k--')
    p.ylim(17,-2)
    p.xlim(110,185)
    p.xticks([125, 150, 175])
    if cnt == 0: p.ylabel('LST [hours]')
    else: p.setp(ax2.get_yticklabels(), visible=False)

    if cnt == 2: p.xlabel('Frequency [MHz]')

fig.subplots_adjust(left=.08, top=.95, bottom=.10, wspace=.3, hspace=.1, right=0.85)
cbar_ax1 = fig.add_axes([0.90, 0.60, 0.03, 0.3])
fig.colorbar(im1, cax=cbar_ax1)
p.xlabel(r'${\rm log}_{10}({\rm Jy})$', fontsize=14)
cbar_ax2 = fig.add_axes([0.90, 0.15, 0.03, 0.33])
fig.colorbar(im2, cax=cbar_ax2)
p.xlabel(r'${\rm radians}$', fontsize=14)

p.show()
