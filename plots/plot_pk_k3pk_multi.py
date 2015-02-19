#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, re, os, glob

o=optparse.OptionParser()
o.add_option('--flux', action='store_true', 
                help='Scale data due to flux calibration errors from pictor. Scales by factor in f option.')
o.add_option('-f', dest='flux_factor', action='store', type='float', default=.736,
            help='scaling factor for flux')
o.add_option('--beam', action='store_true',
                help='scale data by beam square instead of beam')
o.add_option('--afrf', action='store_true',
                help='scale data by factor from aggresive fringe rate filtering.')
o.add_option('-a', dest='afrf_factor', action='store', type='float', default=1.9,
            help='scaling factor for aggressive fring rate filtering.')
o.add_option('--cov', action='store_true',
            help='scale factor for signal loss in covariance removal')
o.add_option('--show', action='store_true',
            help='Show the plot')
opts,args = o.parse_args(sys.argv[1:])
print args

#args = glob.glob('data/final_pspecs/v03*/pspec.npz')
args = glob.glob('data/final_pspecs/v05*/pspec.npz')

ERRORS = {}


fig = p.figure(figsize=(12,7.2))
def dual_plot(kpl, pk, err, pkfold=None, errfold=None, umag=16, f0=.164, color='', ls='.', bins=None, verbose=True, offset=0):
    z = C.pspec.f2z(f0)
    kperp = C.pspec.dk_du(z) * umag #k_perp, from baseline length in wavelengths.
    k = n.sqrt(kpl**2 + kperp**2)
    k3 = n.abs(k**3/(2*n.pi**2))

    #print 'k [h Mpc^-1], P(k) [K^2], err (2sigma)'
    #for _k, _pk, _err in zip(k,pk,err):
    #    print '%6.3f, %9.5f, %9.5f' %(_k, _pk.real/1e6, _err/1e6)
    #print '-'*20

    pk = pk.real #power spectrum is real, so discard the imaginary part.
    p.subplot(121)
    #plot pk with the error bars!
    ebp = p.errorbar(kpl+offset, pk, yerr=err, fmt='.', color=color,capsize=0)#, alpha=.5)
    ebp[-1][0].set_linestyle(ls)

    #for folded power spectrum
    p.subplot(122)
    k0 = n.abs(kpl).argmin() #get index of central k bin.
    if pkfold is None: #if pspec is not folded do a weighted average.
        print 'Folding'
        pkfold = pk[k0:].copy() #take positive k include central bin
        errfold = err[k0:].copy()
        pkpos,errpos = pk[k0+1:].copy(), err[k0+1:].copy()#only positive
        pkneg,errneg = pk[k0-1:0:-1].copy(), err[k0-1:0:-1].copy()#only positive
        pkfold[1:] = (pkpos/errpos**2 + pkneg/errneg**2) / (1./errpos**2 + 1./errneg**2)
        errfold[1:] = n.sqrt(1./(1./errpos**2 + 1./errneg**2))
    else:
        print pkfold.imag
        pkfold = pkfold.real
        
    #plot folded power spectra
    ebd = p.errorbar(k[k0:]+offset, k3[k0:]*pkfold, yerr=k3[k0:]*errfold, fmt='.', color=color, capsize=0)
    ebd[-1][0].set_linestyle(ls)
    if not bins is None:
        _kpls, _k3pks, _k3errs = [], [], []
        for (dn, up) in bins:
            KRMS = False
            if dn < 0 and up > 0: KRMS = True
            ksum,kwgt = 0,0
            dsum,dwgt = 0., 0.
            for _kpl,_pk,_err in zip(kpl, pk, err):
                if dn < _kpl and _kpl <= up:
                    dsum += _pk / _err**2
                    dwgt += 1./ _err**2
                    if KRMS: ksum += _kpl**2 # krms
                    else: ksum += _kpl # kavg
                    kwgt += 1.
            if KRMS: kavg = n.sqrt(ksum/kwgt + kpr**2) # krms
            else: kavg = n.sqrt((ksum/kwgt)**2 + kpr**2) # kavg
            if dwgt == 0: continue
            _pk = dsum / dwgt
            #print _pk
            _err = 1. / n.sqrt(dwgt)
            _k3pk = kavg**3/(2*n.pi**2) * _pk
            _k3err = kavg**3/(2*n.pi**2) * _err
            _kpls.append(dn); _kpls.append(.5*(dn+up)); _kpls.append(up)
            _k3pks.append(_k3pk); _k3pks.append(_k3pk); _k3pks.append(_k3pk)
            _k3errs.append(_k3err); _k3errs.append(_k3err); _k3errs.append(_k3err)
        kpl = n.array(_kpls)
        k3pk = n.array(_k3pks)
        k3err = n.array(_k3errs)
    else:
        kpl = kpl
        k3pk = k3*pk
        k3err = k3*err

    #print k^3/2pi power spectrum, unfolded and filded
    ERRORS[color] = []
    if verbose:
        for _k,_k3pk,_k3err in zip(kpl,k3pk/k3,k3err/k3):
            print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)
            ERRORS[color].append(_k3err)
        print '-'*20
        for _k,_k3pk,_k3err in zip(k[k0:],k3[k0:]*pkfold,k3[k0:]*errfold):
            print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)
        print '-'*20

    print "saving pspec_pk_k3pk.npz"
    print "output @ freq = ",f0
    n.savez('pspec_pk_k3pk.npz',kpl=kpl,pk=pk,err=err,k=k[k0:], k3pk=k3[k0:]*pkfold, k3err=k3[k0:]*errfold,freq=f0)

#######################################################################################                 



FG_VS_KPL_NOS = 168.74e6
FG_VS_KPL = { #K^2
    '-0.049':   1.06015e+15,
    ' 0.000':   8.44e+15, 
    ' 0.049':   1.22195e+15, 
    }

RS_VS_KPL = {} #K^2
RS_VS_KPL_FOLD = {} #K^2
powerspec, powerspecwgt = {}, {}
powerspec_fold, powerspecwgt_fold = {}, {}

#black, blue, orange, green, dark grey
colors = [(0.,0.,0.),  (0, 107/255., 164/255.), (1,128/255.,14/255.), (44/255.,160/255.,44/255.), (89/255.,89/255.,89/255.)]
#lss = ['-','--','-.',':']
lss = ['-','-','-','-']
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    RS_VS_KPL[filename] = {}
    RS_VS_KPL_FOLD[filename] = {}
    powerspec[filename] = {}
    powerspecwgt[filename] = {}
    powerspec_fold[filename] = {}
    powerspecwgt_fold[filename] = {}
    kpl,pk,err = f['kpl'], f['pk'], f['err']
    try: pkfold, errfold = f['pk_fold'],f['err_fold']
    except(KeyError): pkfold,errfold = None, None 
    #save power spectrum per filename. Do not combine
    for _kpl, _pk, _err in zip(kpl, pk, err):
        RS_VS_KPL[filename][_kpl] = (_pk, _err)
        powerspec[filename][_kpl] = _pk / _err**2
        powerspecwgt[filename][_kpl] = 1 / _err**2
    k0 = n.abs(kpl).argmin()
    if not pkfold is None:
        for _kpl, _pk, _err in zip(kpl[k0:], pkfold, errfold):
            RS_VS_KPL_FOLD[filename][_kpl] = (_pk, _err)
            powerspec_fold[filename][_kpl] = _pk / _err**2
            powerspecwgt_fold[filename][_kpl] = 1 / _err**2
#note that there is no need to sum up from different files now.

#freq = f['freq']
freq = .151
#for cnt,filename in enumerate(RS_VS_KPL):
for cnt,filename in enumerate(args):
    print 'loading data from file', filename
    ks = RS_VS_KPL[filename].keys(); ks.sort() 
    d = [powerspec[filename][k]/powerspecwgt[filename][k] for k in ks]
    nos = [1 / n.sqrt(powerspecwgt[filename][k]) for k in ks]

    ks_fold = RS_VS_KPL_FOLD[filename].keys(); ks_fold.sort() 
    d_fold = [powerspec_fold[filename][k]/powerspecwgt_fold[filename][k] for k in ks_fold]
    nos_fold = [1/n.sqrt(powerspecwgt_fold[filename][k]) for k in ks_fold]
    kpl = ks
    kpl_fold = ks_fold
    
    if True:
        print 'Adding Foreground'
        for en, k in enumerate(kpl):
            k = '%6.3f' % k
            if not FG_VS_KPL.has_key(k): continue
            d[en] += FG_VS_KPL[k]
            nos[en] = 2*n.sqrt(FG_VS_KPL_NOS*FG_VS_KPL[k])
        for en,k in enumerate(kpl_fold):
            k = '%6.3f'%k
            if not FG_VS_KPL.has_key(k): continue        
            d_fold[en] += FG_VS_KPL[k]
            nos_fold[en] = 2*n.sqrt(FG_VS_KPL_NOS*FG_VS_KPL[k])

    d,kpl,nos = n.array(d, dtype=n.complex), n.array(ks), n.array(nos)
    d_fold,kpl_fold,nos_fold = n.array(d_fold, dtype=n.complex), n.array(ks_fold), n.array(nos_fold)
    
    if False:
        # PSA32 was calibrated to Pictor A @ 160 MHz = 424 Jy
        # To recalibrate to new Pic A, must multiply by square of ratio of fluxes
        # Jacobs et al 2013 says Pic A = 382 @ 150 MHz, index=-0.76, so at 160 MHz, Pic A = 364 Jy
        #f = 0.76 # psa747 calibration of Pic A = 370.6 Jy @ 160 MHz (which includes resolution effects)
        #f = 0.736 # rescale by (364/424)**2 to correct flux scale
        f = opts.flux_factor
        print 'Scaling data and noise by %f for recalibration to PicA from Jacobs et al. 2013 (PSA32 only)' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    if False:
        f = 2.35 # Use power**2 beam, which is a 1.69/0.72=2.35 penalty factor
        print 'Scaling data and noise by %f for correcting cosmo scalar to use power^2 beam' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    #if True: # For aggressive fringe-rate filtering, change beam area
    if True: # For aggressive fringe-rate filtering, change beam area
        if 'v051' in filename: pass
        #f = opts.afrf_factor
        f = 1.39 # ratio of power**2 beams for filtered * unfiltered beams: 0.306 / 0.162
        print 'Scaling data and noise by %f for beam constriction in aggressive fringe-rate filtering' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f         
    if True: # extra penalty for signal loss in covariance diagonalization
        f = 1.02#updated value
        print 'Scaling data and noise by %f for signal loss in covariance diagonalization' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f

    if True: # penalty for using median statistic
        f = 1/n.log(2)#updated value
        print 'Scaling data and noise by %f for signal loss in covariance diagonalization' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    if True:
        f1 = 1.048
        f2 = 1.013
        f  = 1.0015
        print 'Scaling data and noise by %f(first mode outside horizon), %f(second mode outside horizon), and %f(all modes greater)  for using delay filter.'%(f1,f2,f)
        k0 = n.argmin(n.abs(kpl))
        d[k0+2] *= f1; d[k0-2] *= f1
        d[k0+3] *= f2; d[k0-3] *= f2
        d[k0+4:] *= f; d[:k0-3] *= f
        nos[k0+2] *= f1; nos[k0-2] *= f1
        nos[k0+3] *= f2; nos[k0-3] *= f2
        nos[k0+4:] *= f; nos[:k0-3] *= f

        d_fold[2] *= f1; nos_fold *= f1
        d_fold[3] *= f2; nos_fold *= f2
        d_fold[4:] *= f; nos_fold *= f
    if True:
        f=1 + 2e-9 
        print 'Scaling data and noise by %f for signal loss from flux scale calibration.'%f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
 
    
    if d_fold.size==0: d_fold,nos_fold=None,None
    dual_plot(kpl, d, 2*nos, d_fold, 2*nos_fold, color=colors[0], ls=lss[0], bins=None, f0=freq, offset=.005*cnt)
    colors = colors[1:]+colors[:1]
    print colors
    lss = lss[1:] + lss[:1]


for c in ERRORS:
    print c
    print n.array(ERRORS[c])/1e5
#print n.array(ERRORS['y'])/n.array(ERRORS['c'])
print n.round(ks,8)

tau_h = 100 + 15. #in ns
k_h = C.pspec.dk_deta(C.pspec.f2z(.151)) * tau_h
p.subplot(121)
p.vlines(k_h, -1e7, 1e10, linestyles='--', linewidth=1.5)
p.vlines(-k_h, -1e7, 1e10, linestyles='--', linewidth=1.5)
#p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
p.ylabel(r'$P(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$', fontsize='large')
p.ylim(-.6e7,3e7)
p.grid()


p.subplot(122)
#p.plot([.5], [248**2], 'mv', label='GMRT2013')
#else: p.plot([-.5, .5], [248**2, 248**2], 'mv', label='GMRT2013')
#theoretical 64 thermal noise 
p.vlines(k_h, -1e7, 1e10, linestyles='--', linewidth=1.5)
try: 
    theor_64 = n.load('paper_zsa_lstcnt_zsa_mid_0.164.npz')
    theoks = theor_64['ks'][n.where(theor_64['ks'] < .55)]
    theoerrs = theor_64['T_errs'][n.where(theor_64['ks'] < .55)]
    pktheoerrs = theoerrs*2*n.pi**2 / theoks**3
    pspectheoks = n.concatenate([-1*theoks[::-1], theoks])
    pspectheoerrs = n.concatenate([pktheoerrs[::-1], pktheoerrs])#for 2 sigma and aggressive fringe rate filt.
    p.subplot(121)
    p.plot(pspectheoks, pspectheoerrs, 'k')
    p.plot(pspectheoks, 2*pspectheoerrs, 'k--')
    p.subplot(122)
    p.plot(theoks, theoerrs, 'k')
    p.plot(theoks, 2*theoerrs, 'k--')
except:
    pass
#
#
p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$', fontsize='large')
p.ylim(1e0,1e5)
p.xlim(0, 0.6)
#p.suptitle(colors)
p.grid()
       
        

    
p.show()
