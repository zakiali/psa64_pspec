#! /usr/bin/env python
import numpy as n, pylab as p, capo
from matplotlib import gridspec

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

npz1 = n.load('data/inv_cov_res.npz')
npz2 = n.load('data/inv_cov_fg.npz')
lsts = npz1['lsts'] * 12 / n.pi
fqs = npz1['fqs'] * 1e3
npz = {}
for k in npz1:
    if not k.startswith('x'): continue
    npz[k] = npz1[k]
    npz['fg'+k] = npz2[k]

fig = p.figure(figsize=(9,8))
cols = 7
gs = gridspec.GridSpec(3, cols, width_ratios=[1]*(cols-1)+[.2], height_ratios=[2,1,2]) 

#mode = 'real'
#mx1,mx2,mx3 = 10,15,5
drng1,drng2,drng3 = 20,30,10
mode = 'log'
mx1,mx2,mx3 = 3.5,7,1.
drng1,drng2,drng3=4,8,3

for cnt,bl in enumerate(['fgxbl1','fgxbl2','fgxbl3','xbl1','xbl2','xbl3']):
    if bl.startswith('fg'): scale = 1.#300.
    else: scale = 1.
    x = npz[bl]
    C = cov(x)
    U,S,V = n.linalg.svd(C.conj())
    _C = n.einsum('ij,j,jk', V.T, 1./S, U.T)
    _Cx = n.dot(_C, x)
    #p.plot(S); p.show()
    #print a.miriad.bl2ij(bl), k
    print bl
    #ax1 = p.subplot(2,3,cnt+1)
    ax1 = p.subplot(gs[cnt])
    p.setp(ax1.get_xticklabels(), visible=False)
    if cnt == 0: p.ylabel('LST [hours]')
    else: p.setp(ax1.get_yticklabels(), visible=False)
    im1 = capo.arp.waterfall(x.T, mode=mode, mx=mx1*scale, drng=drng1*scale, extent=(fqs[0],fqs[-1],lsts[-1],lsts[0]))
    p.xticks([148,152,156])
    #p.plot([100,200],[0,0],'k--')
    #p.plot([100,200],[9.5,9.5],'k--')
    #p.ylim(17,-2)
    ax2 = p.subplot(gs[cnt+cols])
    p.setp(ax2.get_xticklabels(), visible=False)
    print C.max(), C.min()
    im2 = capo.arp.waterfall(C, mode=mode, mx=mx2*scale**2, drng=drng2*scale**2, extent=(fqs[0],fqs[-1],fqs[-1],fqs[0]))
    p.xticks([148,152,156])
    p.yticks([148,152,156])
    if cnt == 0: p.ylabel('Frequency [MHz]')
    else: p.setp(ax2.get_yticklabels(), visible=False)

    #ax2 = p.subplot(2,3,3+cnt+1)
    ax3 = p.subplot(gs[cnt+2*cols])
    if cnt == 0: p.ylabel('LST [hours]')
    else: p.setp(ax3.get_yticklabels(), visible=False)
    im3 = capo.arp.waterfall(_Cx, mode=mode, mx=mx3, drng=drng3, extent=(fqs[0],fqs[-1],lsts[-1],lsts[0]))
    p.xticks([148,152,156])
    if cnt == 3: p.xlabel('Frequency [MHz]')
    #p.subplot(311); capo.arp.waterfall(x, mode='real', mx=10,drng=20); p.colorbar()
    #p.subplot(323); capo.arp.waterfall(C)
    #p.subplot(324); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
    #p.subplot(313); capo.arp.waterfall(_Cx, mode='real', mx=5, drng=10); p.colorbar()

fig.subplots_adjust(left=.1, top=.95, bottom=.10, wspace=.2, hspace=.1, right=0.9)
#cbar_ax1 = fig.add_axes([0.90, 0.65, 0.03, 0.3])
cbar_ax1 = p.subplot(gs[cols-1])
fig.colorbar(im1, cax=cbar_ax1, ticks=[0,1,2,3], format='%+d')
cbar_ax1.yaxis.set_label_position("right")
p.ylabel(r'${\rm log}_{10}({\rm mK})$', fontsize=13)

#cbar_ax2 = fig.add_axes([0.90, 0.35, 0.03, 0.3])
cbar_ax2 = p.subplot(gs[2*cols-1])
fig.colorbar(im2, cax=cbar_ax2, ticks=[0,2,4,6], format='%+d')
cbar_ax2.yaxis.set_label_position("right")
p.ylabel(r'${\rm log}_{10}({\rm mK}^2)$', fontsize=13)

#cbar_ax3 = fig.add_axes([0.90, 0.10, 0.03, 0.33])
cbar_ax3 = p.subplot(gs[3*cols-1])
fig.colorbar(im3, cax=cbar_ax3, ticks=[-2,-1,0,1], format='%+d')
cbar_ax3.yaxis.set_label_position("right")
p.ylabel(r'${\rm log}_{10}({\rm mK}^{-1})$', fontsize=13)
p.show()
