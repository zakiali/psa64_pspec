#! /usr/bin/env python
import numpy as n, pylab as p

def Tcmb(z):
    return 2.7e3 * (1+z)

def Tgas(z):
    return Tcmb(z) * ((1+z)/150.)

def T0(z):
    return 26.7*((1+z)/10)**.5

def k3pk_patchy(T_b, x_i=0.5, kmax_kmin=2):
    x_H = 1-x_i
    return (x_H - x_H**2) * T_b**2 / n.log(kmax_kmin)

def Ts(z,tb):
    return Tcmb(z)*(1-(tb/T0(z)))**-1

def Ts_post_2s(deli,pks,errs,z):
    #given a range of ts, find the posterior dist and return 2 sigma upper limit.   
    tss = n.linspace(1000, 20000, 5000) #range of ts
    #assuming tbs are positive
    tbs = -1*T0(z)*(1 - (Tcmb(z)/tss))
    s = tbs**2 * deli
    #get posterior for del2 and compare to tbs to map back to tss 2 sigma upper limit. 
    data = []
#    s = n.linspace(-50, 1000, 1000)
    for ss in s:
        data.append(n.exp(-.5*n.sum((pks.real-ss)**2 / errs**2)))
    data = n.array(data)
#    p.plot(s,data)
#    p.show()
    data /= n.max(data)
    data_c  = n.cumsum(data)
    data_c /= data_c[-1]
    s2hi = s[data_c>1-.0227][0]

    tbs_f = n.sqrt(s2hi/deli)
    #print tbs_f
    #print tbs
    ts_2sig = tss[tbs<tbs_f][0]
    print 'tspin', ts_2sig #tss is a decreasing function
    
#    print tbs 
#    print 'tbs_f ', tbs_f 
    
    return ts_2sig
        

limits = {}
#limits['.05-.1'] = (n.arange(.05, .1, .001), 3100.,'c')
#limits['.1-.2'] = (n.arange(.1, .2, .001), 1450.,'m')
#limits['.2-.4'] = (n.arange(.2, .4, .001), 3460.,'b')
##limits['.4-.5'] = (n.arange(.4, .52, .001), 8210., 'y')

# Errors are 2sigma
#if False: # ver1 limits
#    limits['0.081'] = (0.081,11757, 1791, 'k')
#    limits['0.108'] = (0.108, 1943,  748, 'k')
#    limits['0.134'] = (0.134, 2330, 1256, 'k')
#    limits['0.161'] = (0.161, 2818, 1099, 'k')
#    limits['0.188'] = (0.188, 2920, 1902, 'k')
#    limits['0.215'] = (0.215, 3256, 2207, 'k')
#    limits['0.241'] = (0.241, 1882, 3219, 'k')
#    limits['0.268'] = (0.268, 2161, 4670, 'k')
#    limits['0.295'] = (0.295, 5741, 6295, 'k')
#    limits['0.322'] = (0.322, 3727, 7530, 'k')
#    limits['0.349'] = (0.349, 3690, 9600, 'k')
#    limits['0.376'] = (0.376, 9466,14191, 'k')
#    limits['0.402'] = (0.402,10451,18724, 'k')
#    limits['0.429'] = (0.429, 6751,19389, 'k')
#    limits['0.456'] = (0.456, 5531,22167, 'k')
#    limits['0.483'] = (0.483,26606,30157, 'k')
#    limits['0.510'] = (0.510,61758,35445, 'k')
#elif False: # v002 limits, from v009
#    limits['0.082'] = (0.082, 5480,  911, 'k')
#    limits['0.108'] = (0.108, 4085, 1641, 'k')
#    limits['0.135'] = (0.135, 3373, 2038, 'k')
#    limits['0.162'] = (0.162, 2970, 2496, 'k')
#    limits['0.188'] = (0.188, 2158, 3465, 'k')
#    limits['0.215'] = (0.215, -127, 4853, 'k')
#    limits['0.242'] = (0.242, -582, 6336, 'k')
#    limits['0.269'] = (0.269,-5041, 7488, 'k')
#    limits['0.295'] = (0.295,-4808,10178, 'k')
#    limits['0.322'] = (0.322, 9408,14388, 'k')
#    limits['0.349'] = (0.349,19362,18826, 'k')
#    limits['0.376'] = (0.376, 3561,22972, 'k')
#    limits['0.403'] = (0.403,-8184,26199, 'k')
#    limits['0.429'] = (0.429,-4183,31652, 'k')
#    limits['0.456'] = (0.456,  200,37778, 'k')
#    limits['0.483'] = (0.483,-5975,44836, 'k')
#    limits['0.510'] = (0.510,33556,61554, 'k')
#else: # v006 limits, bl28_pg4_iter9_autox
#    limits['0.082'] = (0.082, 6457,   670, 'k')
#    limits['0.108'] = (0.108, 5999,  2257, 'k')
#    limits['0.135'] = (0.135,14789,  6404, 'k')
#    limits['0.162'] = (0.162,11582,  3526, 'k')
#    limits['0.188'] = (0.188, 3852,  2322, 'k')
#    limits['0.215'] = (0.215, 3323,  1625, 'k')
#    limits['0.242'] = (0.242, 2448,  1677, 'k')
#    limits['0.269'] = (0.269, -506,  2164, 'k')
#    limits['0.295'] = (0.295, 1003,  2333, 'k')
#    limits['0.322'] = (0.322,11498,  3406, 'k')
#    limits['0.349'] = (0.349,13086,  5855, 'k')
#    limits['0.376'] = (0.376, 8773,  7753, 'k')
#    limits['0.403'] = (0.403, 2049,  9535, 'k')
#    limits['0.429'] = (0.429,-3360,  7591, 'k')
#    limits['0.456'] = (0.456,  -75, 10571, 'k')
#    limits['0.483'] = (0.483, 6539, 10439, 'k')
#    limits['0.510'] = (0.510,83498, 26917, 'k')

#PSA64 limits. 
if True:
    limits['.100'] = (0.100, 435.60, 102.88, 'k')
    limits['.148'] = (0.148, 863.10, 195.78, 'k')
    limits['.197'] = (0.197, 178.14, 203.66, 'k')
    limits['.246'] = (0.246, 266.52, 376.30, 'k')
    limits['.295'] = (0.295, 303.51, 551.55, 'k')
    limits['.345'] = (0.345, 1177.44, 831.86, 'k')
    limits['.394'] = (0.394, -1288.74, 1265.40, 'k')
    limits['.443'] = (0.443, -2.25, 1487.52, 'k')
    limits['.492'] = (0.492, 2083.53, 2562.52, 'k')

d2points = []
d2errs = []
ks = []
for L in limits:
    ks.append(limits[L][0])
    d2points.append(limits[L][1])
    d2errs.append(limits[L][2])
ks = n.array(ks)
d2points = n.array(d2points)
d2errs = n.array(d2errs)

axes = []
fig = p.figure(figsize=(10,4))
#kwids = n.arange(1.1,3000.,.1)
#kwids = 10**n.arange(0.1, 6, .1)
#kcens = 10**n.arange(-2, 2, .1)
kmins = 10**n.arange(-3, 2.1, .03)
kmaxs = 10**n.arange(-3, 2.1, .03)
#for xi,sty in [(.5,'-'), (.3,'--'), (.1,'-.'), (.03,':')]:
for cnt,(xi,sty) in enumerate([(.1,'-'), (.3,'--'), (.5,'-.')]):
    T_b_lim = []
    xs = [kmins] * len(kmaxs); xs = n.array(xs).transpose()
    ys = [kmaxs] * len(kmins); ys = n.array(ys)
#    print xs.shape, ys.shape
    for kmax in kmaxs:
        T_b_lim_kwid = []
        for kmin in kmins:
            tmp = []
            include = n.logical_and(ks>kmin, ks<kmax)
            if n.all(n.logical_not(include)): 
                T_b_lim_kwid.append(1.576) #K
                continue
#            print kmin, kmax, ks
#            print include
            lim = d2points[include]
#            print lim
            err = d2errs[include]
            avg = n.average(k3pk_patchy(1., xi, kmax/kmin))
            sig2 = Ts_post_2s(avg, lim, err,8.4) 
            T_b_lim_kwid.append(sig2)
            
                   
 
#            for L in limits:
#                ks, k3pk, err, clr = limits[L]
#                #if ks < kcen - kwid or ks > kcen + kwid:
#                if ks < kmin or ks > kmax:
#                    tmp.append(391)
#                    continue
#                lim = k3pk + err
#                ks = n.linspace(ks-.014, ks+.14, 100)
#                #avg = n.average(k3pk_patchy(ks, 1., xi, kwid))
#                avg = n.average(k3pk_patchy(1., xi, kmax/kmin))
#                tmp.append(n.sqrt(lim / avg))
#            T_b_lim_kwid.append(min(tmp+[391]))
        T_b_lim.append(T_b_lim_kwid)
    T_b_lim = n.array(T_b_lim)
#    T_s_lim = Ts(8.4, -1*T_b_lim)
    print T_b_lim.shape
#    print T_s_lim.shape
    #p.plot(Rs, T_b_lim, 'k'+sty)
    #p.fill_between(kwids, T_b_lim, 1000*n.ones_like(kwids), facecolor='k', alpha=0.2)
#    p.subplot(1,3,cnt+1)
    axes.append(p.subplot(1,3,cnt+1))
    #image = axes[cnt].imshow(T_b_lim, vmax=150, vmin=0, interpolation='nearest', aspect='auto', origin='lower',
    image = axes[cnt].imshow(T_b_lim/1000., vmax=10, vmin=0, interpolation='nearest', aspect='auto', origin='lower',
        #extent=(-3,2,-3,2))
        extent=(1e-3,1e2,1e-3,1e2))
    p.xscale('log')
    p.yscale('log')
    p.xlim(1e-3,0.45)
    p.ylim(1e-1,1e2)
    #p.fill_between([-3,3],[-3,3],[-3,-3], facecolor='w', alpha=1)
    if xi == 0.5: p.title('$x_i=%3.1f$'%xi)
    else: p.title('$x_i=%3.1f,%3.1f$' % (xi,1-xi))
    #p.contourf(ys, xs, T_b_lim, [50,100,150,200,250,300,350,400])
    p.fill_between([1e-3,1e3],[1e-3,1e3],[1e-3,1e-3], facecolor='w', alpha=1)
    #p.colorbar()
    p.xlabel(r'$k_{\rm min}$ [$h^{-1}\ {\rm Mpc}$]', size=14)
    if cnt == 0: p.ylabel(r'$k_{\rm max}$ [$h^{-1}\ {\rm Mpc}$]', size=14)
    p.grid()
#for the circle and triangle
#    if cnt in [1,2]:
#        p.plot([.2], [30.], 'ko')
#        p.plot([.1], [30.], 'k^')
        #p.plot([.3], [10.], 'ko')
        #p.plot([.15], [10.], 'k^')
    #if cnt == 2: p.colorbar()
    
#p.colorbar()
cbar_axis=fig.add_axes([.92,.15,.03,.75])
fig.colorbar(image, cax=cbar_axis)
p.subplots_adjust(left=.1,bottom=.15,top=.90,right=.90)
p.show()

exit()#XXX
#p.semilogx([1., 1000], [400, 400], 'k--')
p.fill_between([1., 1e8], [400, 400], [1000,1000], facecolor='m', alpha=.5)
p.semilogx([1., 1e8], [30, 30], 'k--')
p.ylim(0,500)
p.xlim(1e1,1e6)
p.ylabel(r'$|\langle T_b\rangle|$ [${\rm mK}$]', size=14)
p.xlabel(r'$k_{\rm max}/k_{\rm min}$ [$h^{-1}\ {\rm Mpc}$]', size=14)
p.grid()

for L in limits:
    ks, k3pk, err, clr = limits[L]
    lim = k3pk + err
    kmax,kmin = 10.,.1
    for xi in [.3, .5]:
        avg = n.average(k3pk_patchy(1., xi, kmax/kmin))
        print L, xi, n.sqrt(lim / avg)

if True: # set true to get the scalings for the lidz curves
    import glob, scipy.interpolate, re, os, capo as C, scipy.special

    re_z = re.compile(r'power_21cm.*_z(\d+\.\d+).*\.dat')
    # xi = .54, .71
    #files_norm = [glob.glob('lidz_mcquinn_k3pk/power*%s*.dat' % s)[0] for s in ['7.3','7.0']]
    files_norm = [glob.glob('lidz_mcquinn_k3pk/power*%s*.dat' % s)[0] for s in ['8.34', '8.76']]
    # xi = .77, .63, .51, .41
    # xi = .51, .77
    #files_hmass = [glob.glob('lidz_mcquinn_k3pk/hmass/power*%s*.dat' % s)[0] for s in ['7.1','7.3','7.4','7.6']]
    #files_hmass = [glob.glob('lidz_mcquinn_k3pk/hmass/power*%s*.dat' % s)[0] for s in ['7.4','7.1']]
    files_hmass = [glob.glob('lidz_mcquinn_k3pk/hmass/power*%s*.dat' % s)[0] for s in ['8.34','8.76']]

    #for cnt,filename in enumerate(glob.glob('lidz_mcquinn_k3pk/*dat')):
    for cnt,files in enumerate([files_norm, files_hmass]): #glob.glob('lidz_mcquinn_k3pk/hmass/power_21cm_hmass_*dat')):
        clr = 'cm'[cnt]
        for cnt,filename in enumerate(files):
            sym = ['-','--','-.',':'][cnt]
            print 'Reading', filename, clr + sym
            d = n.array([map(float, L.split()) for L in open(filename).readlines()])
            ks, pk = d[:,0], d[:,1]
            z_file = float(re_z.match(os.path.basename(filename)).groups()[0])
            z = C.pspec.f2z(.160)
            k3pk = ks**3 / (2*n.pi**2) * pk
            mdl = scipy.interpolate.interp1d(ks, k3pk, kind='linear')
            Ls = limits.keys(); Ls.sort()
            ks = n.array([limits[L][0] for L in Ls])
            k3pk = n.array([limits[L][1] for L in Ls])
            err = n.array([limits[L][2] for L in Ls])
            lvl = mdl(ks)
            Ts = n.sqrt((k3pk+err)/lvl)
            print ks
            print n.around(Ts)
            Tmin = Ts.min()
            print '-'*10
            p.plot([1., 1000], [Tmin, Tmin], clr+sym)

p.show()
import sys; sys.exit(0)

def V_to_pk(V, L, nbins=200):
    dL = L / V.shape[0]
    kmin,kmax = 2*n.pi/L, 2*n.pi/dL
    P_k = n.abs(n.fft.ifftn(V)/dL**3)**2 * L**3
    print P_k[0,0,0], L, dL
    kx = 2*n.pi * n.fft.fftfreq(V.shape[0], d=dL); kx.shape = (kx.size,1,1)
    ky = 2*n.pi * n.fft.fftfreq(V.shape[1], d=dL); ky.shape = (1,ky.size,1)
    kz = 2*n.pi * n.fft.fftfreq(V.shape[2], d=dL); kz.shape = (1,1,kx.size)
    k = n.sqrt(kx**2 + ky**2 + kz**2)
    print k[0,0,0], k[1,0,0]
    hist_sum,bins = n.histogram(k, range=(kmin,kmax), bins=nbins, weights=P_k)
    hist_wgt,bins = n.histogram(k, range=(kmin,kmax), bins=nbins, weights=n.ones_like(P_k))
    hist = hist_sum / n.where(hist_wgt > 0, hist_wgt, 1)
    return 0.5*(bins[:-1] + bins[1:]), hist
    
kmin,kmax = .03, 3. # h Mpc^-1
L,dL = 2*n.pi/kmin, 2*n.pi / kmax
SZ = int(L / dL)

print L, dL
print SZ


# Add in bubbles
def ionize(V, L, xi=0.5, avg_log10R=n.log10(2*n.pi/0.25), sig_log10R=.1):
    x,y,z = n.indices(V.shape)
    print n.average(V)
    while n.average(V) > xi:
        r = 10**n.random.normal(loc=avg_log10R, scale=sig_log10R)
        r /= L / V.shape[0]
        print n.average(V), r
        x0,y0,z0 = n.random.randint(V.shape[0], size=len(V.shape))
        print x0, y0, z0
        V = n.where((x-x0)**2+(y-y0)**2+(z-z0)**2 <= r**2, 0, V)
    print n.average(V)
    return V

def ionize_square(V, L, xL):
    x,y,z = n.indices(V.shape, dtype=n.float32)
    xL = xL / (L / V.shape[0])
    x = n.where((x/xL) % 2 < 1, 0, 1)
    y = n.where((y/xL) % 2 < 1, 0, 1)
    z = n.where((z/xL) % 2 < 1, 0, 1)
    xyz = (x + y + z) % 2
    return V * xyz
        
#for xL in n.arange(n.log10(kmin), n.log10(kmax), .1):
#    V = n.ones((SZ,SZ,SZ), dtype=n.float32)
#    xL = 2*n.pi/(10**xL)
#    print xL, dL
#    V = ionize_square(V, L, xL)
#    print n.average(V)
#    ks, pk = V_to_pk(V, L)
#    p.plot(ks, ks**3/(2*n.pi**2) * pk)
#    print 'done'
V = n.ones((SZ,SZ,SZ), dtype=n.float32)
V = n.random.normal(scale=V)
V *= (Tcmb(8.) - Tgas(8.))
ks, pk = V_to_pk(V, L)
p.plot(ks, ks**3/(2*n.pi**2) * pk)

p.show()

