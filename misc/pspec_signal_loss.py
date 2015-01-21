import numpy as n
import pylab as p
import aipy as a
import optparse, sys, capo, glob

o = optparse.OptionParser()
a.scripting.add_standard_options(o,chan=True)
o.add_option('--pk', type='float', default=300.6,
    help='level of flat P(k).')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV('../plots/data/binned_uv/lst.2456242.51881.uvA')
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])

afreqs = freqs.take(chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

jy2T = capo.pspec.jy2T(fq)

B = sdf * afreqs.size
etas = n.fft.fftshift(capo.pspec.f2eta(afreqs))
kpl = etas * capo.pspec.dk_deta(z)

bm = n.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq) * 2.35
scalar = capo.pspec.X2Y(z) * bm * B
print 'Freq:',fq
print 'z:', z
print 'B:', B
print 'scalar:', scalar

def noise(nsamp, var):
    return n.random.normal(size=nsamp)*n.sqrt(var/2.) + 1j*n.random.normal(size=nsamp)*n.sqrt(var/2. )

if False:
    #This is to figure out what amplitude noise we need in the input to the
    #power spectrum to getthe output power value. Note that this is a one sample
    #number. Divide by root N to get actual output.  
    pk = opts.pk
    variance_in_tau = pk * bm*B / scalar / jy2T**2
    print 'variance of tau',variance_in_tau

    data_in_tau = noise(203, variance_in_tau)
    #print data_in_tau
    data_in_freq = n.fft.ifft(data_in_tau)

    print 'var in tau,', n.average(n.abs(data_in_tau)**2) #- n.average(n.abs(data_in_tau))**2
    print 'var in freq,', n.var(data_in_freq)
    normalization = n.var(data_in_freq)/n.var(data_in_tau)
    print 'normaliztion', normalization


#This is actual signal loss. Compares covariance weighted data (which has signal
#loss) and compare it with non covariance weighted data.
if True:
    boots = {}
    for dir in args:
        if not dir in boots: 
            boots[dir] = {}
            boots[dir]['pk_vs_t'] = {}
            boots[dir]['nocov_vs_t'] = {}
        for bfile in glob.glob(dir+'/pspec_boot00*'):
            bootnumber = bfile.split('/')[-1]
            print dir+'/'+bootnumber
            boot_npz = n.load(bfile)
            pk_vs_t = boot_npz['pk_vs_t']
            nocov_vs_t = boot_npz['nocov_vs_t']
            boots[dir]['pk_vs_t'][bootnumber] = pk_vs_t
            boots[dir]['nocov_vs_t'][bootnumber] = nocov_vs_t

    print boots.keys()


    nboots = len(boots[dir]['pk_vs_t'])
    print nboots
    pk_vs_t_avg = {}
    nocov_vs_t_avg = {}
    signal_loss = {}
    bns = boots[dir]['pk_vs_t'].keys()
    print bns
    for key in boots.keys():
        if not key in pk_vs_t_avg: 
            pk_vs_t_avg[key] = n.zeros(boots[key]['pk_vs_t'][bns[0]].shape, dtype = n.complex)
            nocov_vs_t_avg[key] = n.zeros(boots[key]['nocov_vs_t'][bns[0]].shape, dtype = n.complex)
        for bootnumber in bns:
            pk_vs_t_avg[key] += boots[key]['pk_vs_t'][bootnumber]
            nocov_vs_t_avg[key] += boots[key]['nocov_vs_t'][bootnumber]

        pk_vs_t_avg[key]/=20.
        nocov_vs_t_avg[key]/=20.
#        sl = nocov_vs_t_avg[key] / pk_vs_t_avg[key]
#        signal_loss[key] = n.median(sl, axis=-1)
        #signal_loss[key] = n.median(nocov_vs_t_avg[key], axis=-1) / n.median(pk_vs_t_avg[key], axis=-1)
        signal_loss[key] = n.median(nocov_vs_t_avg[key], axis=-1) - n.median(pk_vs_t_avg[key], axis=-1)

    print signal_loss[key]
#    print pk_vs_t_avg
#    print nocov_vs_t_avg

#    p.figure(1)
#    capo.arp.waterfall(pk_vs_t_avg[key]);p.colorbar(shrink=.5)
#    p.figure(2)
#    capo.arp.waterfall(nocov_vs_t_avg[key]);p.colorbar(shrink=.5)

    kpl = boot_npz['kpl']

    for key in boots.keys():
        print key
#        p.figure(3)
#        p.plot(signal_loss[key].imag)
        p.figure(4)
        p.plot(kpl, signal_loss[key].real, 'o', label=key)
        

    p.legend()
    p.show()
