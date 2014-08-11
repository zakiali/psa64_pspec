#! /usr/bin/env python
'''Uses lst binned data with weights'''
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse

args = sys.argv[1:]
args = ['data/raw/omnical_lstbin_v2Bf_I.npz','data/raw/omnical_lstbin_v2_I.npz']

sep = 'sep6'
lst_rng = (3,5)


jy2T = {}
jy2T_df = {}
dat = {}; wgt = {}
raw_dat, raw_wgt = {}, {}
df_dat, df_wgt = {}, {}
dt_dat, dt_wgt = {}, {}
#raw_dat, raw_wgt = 0., 0.
#df_dat, df_wgt = 0., 0.
#dt_dat, dt_wgt = 0., 0.
for filename in args:
    print 'Reading', filename
    ftype = filename.split('.')[-2]
    npz = n.load(filename)
    lsth = npz['lsts']*12/n.pi
    lstrng = n.logical_and(lsth>lst_rng[0], lsth<lst_rng[-1])
    print 'B =', .1/203
    if not jy2T.has_key(ftype):
        freqs = n.arange(.1,.2,.1/203) 
        freqs_df = 0.5 * (freqs[1:] + freqs[:-1])
        jy2T[ftype] = C.pspec.jy2T(freqs)
        jy2T_df[ftype] = C.pspec.jy2T(freqs_df)
    #get data and wgt
    dat = npz[sep][lstrng,:]#*npz['wgt_'+sep]
    bwgt = dat != 0.
    bwgt = bwgt.astype(int)
    
    # Re-calibrate to new picA flux scale.XXX need to do this but to something else.
    #d *= (363.6/424) * (freqs/.160)**(-.76+.95)
    d_dt = (dat[1:,] - dat[:-1,]) * n.logical_and(bwgt[1:,], bwgt[:-1,])
    w_dt = 2*n.logical_and(bwgt[1:,], bwgt[:-1,])
    d_df = dat[:,1:] - dat[:,:-1] * n.logical_and(bwgt[:,1:], bwgt[:,:-1])
    w_df = 2*n.logical_and(bwgt[:,1:], bwgt[:,:-1])

    #Average all of these in time. thats the sum. rms os the square
    raw_dat[ftype] = raw_dat.get(ftype,0) + n.sum(n.abs(dat)**2, axis=0)
    raw_wgt[ftype] = raw_dat.get(ftype,0) + n.sum(bwgt, axis=0)
    df_dat[ftype] = df_dat.get(ftype,0) +   n.sum(n.abs(d_df)**2, axis=0)
    df_wgt[ftype] = df_wgt.get(ftype,0) +   n.sum(w_df, axis=0)
    dt_dat[ftype] = dt_dat.get(ftype,0) + n.sum(n.abs(d_dt)**2, axis=0)
    dt_wgt[ftype] = dt_wgt.get(ftype,0) + n.sum(w_dt, axis=0)
    #raw_dat[ftype] = raw_dat.get(ftype,0) + n.abs(n.sum(dat, axis=0))**2
    #raw_wgt[ftype] = raw_dat.get(ftype,0) + n.sum(bwgt, axis=0)
    #df_dat[ftype] = df_dat.get(ftype,0) +   n.abs(n.sum(d_df, axis=0))**2
    #df_wgt[ftype] = df_wgt.get(ftype,0) +   n.sum(w_df, axis=0)
    #dt_dat[ftype] = dt_dat.get(ftype,0) + n.abs(n.sum(d_dt, axis=0))**2
    #dt_wgt[ftype] = dt_wgt.get(ftype,0) + n.sum(w_dt, axis=0)

    print 't'


for i,ftype in enumerate(jy2T.keys()):
    raw = n.sqrt(n.where(raw_wgt[ftype] > 0, raw_dat[ftype] / raw_wgt[ftype], 0)) * jy2T[ftype]
    df = n.sqrt(n.where(df_wgt[ftype] > 0, df_dat[ftype] / df_wgt[ftype], 0)) * jy2T_df[ftype]
    dt = n.sqrt(n.where(dt_wgt[ftype] > 0, dt_dat[ftype] / dt_wgt[ftype], 0)) * jy2T[ftype]
    p.semilogy(freqs, raw, 'k')
    p.semilogy(freqs_df, df, 'm')
    p.semilogy(freqs, dt, 'c')
p.grid()
p.xlim(.115,.187)
p.ylim(1,3e3)
p.xlabel('Frequency [GHz]')
p.ylabel('Brightness Temperature [mK]')
p.savefig('noise_t_35.png', format='png')
p.show()
