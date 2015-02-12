import numpy as n, pylab as p, aipy as a
import capo as C

aa = a.cal.get_aa('psa6240_v003', 0.0004926108374384237, .1, 203)
bls, conj = C.red.group_redundant_bls(aa.ant_layout)
#dont use north south baselines in beamform.
dont_use = ['%d,0'%i for i in range(1,8)]
ubls = bls.keys()
for i in dont_use: ubls.remove(i)
print len(ubls)

Gianni=True

POLY = [-551229901929.10315, 685757303983.00293, -370738356723.81897, 113784843213.44109, -21686535302.711388, 2628544446.213274, -197867884.08778715, 8457793.3780113906, -157170.17280964542]
#Gianni's bandpass
if Gianni:
    POLY = [5.65900306e-14,  -7.20008736e-11,   4.05864355e-08, -1.33028779e-05,   2.79385101e-03,  -3.89867874e-01, 3.61453028e+01,  -2.14668637e+03,   7.40998469e+04, -1.13249078e+06]
NCHAN = 203
CNT=3000
window = a.dsp.gen_window(NCHAN, 'blackman-harris')
fq = aa.get_afreqs()
delays = n.fft.fftfreq(NCHAN, fq[1]-fq[0])
dk_deta = C.pspec.dk_deta(C.pspec.f2z(.1515))
BP = n.polyval(POLY, fq)
if Gianni:
    BP = 1/BP
offset = 1e5
nmodes = 1017

avg1 = 0
avg2 = 0
avg = 0
for i in xrange(CNT):
    ns = []
    #for nbl in range(len(ubls)):
    for nbl in range(nmodes):
        #create new random noise for a new type of baseline.
        ns.append(n.random.normal(size=NCHAN) * n.exp(2j*n.pi*n.random.uniform(size=NCHAN)))
#    ns = n.random.normal(size=NCHAN) * n.exp(2j*n.pi*n.random.uniform(size=NCHAN))
    ns = n.array(ns)
    ns1 = BP *  ( offset + ns.copy())
    ns2 = BP *  ( offset + ns.copy())
    ns1_avg = n.average(ns1, axis = 0)
    ns2_avg = n.average(ns2, axis = 0)
    poly = n.polyfit(fq, ns2_avg, deg=9)
    bp_eval = n.polyval(poly, fq)
    ns1 = offset * (ns1 / (BP * offset) - 1)
    ns2 = offset * (ns2 / bp_eval - 1)
#    print n.abs(ns2-ns1)
    
    _ns1 = n.fft.fft(window*ns1)
    _ns2 = n.fft.fft(window*ns2)
#    p.plot(ns1,'k')
#    p.plot(ns2,'b')
#    p.show()
    avg1 += n.abs(_ns1)**2
    avg2 += n.abs(_ns2)**2
    avg += n.abs(_ns2)**2/n.abs(_ns1)**2

avg1/=CNT
avg2/=CNT
avg/=CNT

p.semilogy(delays*dk_deta, n.abs(n.average(avg,axis=0)-1), '.')
p.plot(n.array([-115,-115])*dk_deta,[1e-10,1e10], 'k--')
p.plot(n.array([ 115, 115])*dk_deta,[1e-10,1e10], 'k--')
p.xlim(-.5,.5)
p.ylim(1e-5,1e1)
p.show()
