#! /usr/bin/env python
import capo as C
import aipy as a, numpy as n, pylab as p

def get_data(filenames, antstr, polstr, verbose=False):
    # XXX could have this only pull channels of interest to save memory
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
            #if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            #pol = a.miriad.pol2str[uv['pol']]
            #if not dat[bl].has_key(pol):
            #    dat[bl][pol],flg[bl][pol] = [],[]
            #dat[bl][pol].append(d)
            #flg[bl][pol].append(f)
    return n.array(lsts), dat, flg



aa = a.cal.get_aa('psa6240_v003', .1/203, .1, 203)
chans = n.array([120])
freqs = aa.get_afreqs()
freqs = freqs.take(chans)
import glob
times, data, flags = get_data(glob.glob('data/lstdata/uvA_fg/lst*uvAA'), 'cross', 'I')
vis  = n.array(data[283])
wgts = n.array(flags[283])
vis  = vis.take(chans, axis=1)
wgts = n.logical_not(wgts.take(chans, axis=1))
print vis.shape, wgts.shape

#unfringe rate filtered visibility
tw = a.dsp.gen_window(vis.shape[0], window='blackman-harris'); tw.shape += (1,)
vis_tft = n.fft.fftshift(n.fft.ifft(vis*tw,axis=0),axes=0)


frates = n.fft.fftshift(n.fft.fftfreq(len(vis), 42.8) * 1e3) #mHz
step = frates[1] - frates[0]
f1, f2 = frates[0] - .5*step, frates[-1]+.5*step
fq = n.average(freqs)

beam_w_fr = C.frf_conv.get_beam_w_fr(aa,(0,26))
tbins, firs, frbins, frspace = C.frf_conv.get_fringe_rate_kernels(beam_w_fr, 42.8, 401)
fr_bins = n.fft.fftshift(n.fft.fftfreq(401, 42.8))

bwfr10 = beam_w_fr.take(chans, axis=0)[0](fr_bins)
firs = firs.take(chans, axis=0).T
frspace= frspace.take(chans, axis=0).T
firs = firs.conj()

vis_frf = n.zeros_like(vis)
wgts_frf = n.zeros_like(vis)
filtered = n.zeros_like(vis)
for ch in xrange(len(chans)):
    #vis_frf[:,ch] = n.convolve(wgts[:,ch]*vis[:,ch], firs, mode='same')
    #wgts_frf[:,ch] = n.convolve(wgts[:,ch], n.abs(firs), mode='same')
    #filtered[:,ch] = n.where(wgts_frf[:,ch]>0, vis_frf[:,ch]/wgts_frf[:,ch], 1)
    vis_frf[:,ch] = n.convolve(wgts[:,ch]*vis[:,ch], firs[:,ch], mode='same')
    wgts_frf[:,ch] = n.convolve(wgts[:,ch], n.abs(firs[:,ch]), mode='same')
    filtered[:,ch] = n.where(wgts[:,ch]>0, vis_frf[:,ch]/wgts_frf[:,ch], 1)
ww = n.where(wgts_frf>.1, 1, 0)
filtered = ww*filtered

filtered_in_fr = n.fft.fftshift(n.fft.ifft(filtered*tw,axis=0)) #tw is a windowing function.

frbins = frbins*1000

rat =  n.max(n.abs(vis_tft))/n.max(n.abs(filtered_in_fr))
print rat
resc_index = n.where( n.abs((n.abs(frspace) - (1/rat))) < .03 )[0][0]
rescale = n.max(n.abs(vis_tft))/rat
#frat = n.max(n.abs(vis_tft)) * rat
print n.abs(frspace[resc_index])

p.subplot(111)
p.plot(frates, n.abs(vis_tft), label='Unfiltered')
p.plot(frbins, n.abs(frspace)*rescale*rat, label='Filter')
p.plot(frates, n.abs(filtered_in_fr), label='Filtered')
p.xlim(-.6,1.5)
p.xlabel('Fringe Rate (milli Hz)', size='large')
p.ylabel('Amplitude (Jy)', size='large')
p.legend(loc=2)


p.show()
