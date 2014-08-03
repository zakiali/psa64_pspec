import frf_conv
import aipy as a, numpy as n, pylab as p
import capo as C

#get antenna array for frf filters.
aa = a.cal.get_aa('psa6240_v003', .1/203, .1, 203)
#channel range
chans = n.array([120])
freqs = aa.get_afreqs()
freqs = freqs.take(chans)

#load data
data = n.load('data/raw/omnical_lstbin_v2_I.npz')
#data = n.load('data/raw/omnical_lstbin_v2Bf_I_0.npz')
#data = n.load('data/raw/omnical_lstbin_v2_xx.npz')
sep = 'sep24' #the 0,1 spacing 

#get vis, and wgts data in the given channel range and put into observation time order.
vis = data[sep]
wgts = data['wgt_'+sep] / 1e-6
lsts = data['lsts']
vis = vis.take(chans, axis=1)
wgts = wgts.take(chans, axis=1)
vis = n.concatenate([vis[-200:], vis[:1460]], axis=0)#[300:800] #want in time order for fr filter
wgts = n.concatenate([wgts[-200:], wgts[:1460]], axis=0)#[300:800] #want in time order for frf
lsts = n.concatenate([lsts[-200:], lsts[:1460]], axis=0)#[300:800] #want in time order for frf

#for plotting the correct lsts.
step = lsts[1] - lsts[0]
if lsts[0] > lsts[-1]:#wrapped
    diff = 2*n.pi - lsts[0]
    lstsmod = ((lsts + diff)%(2*n.pi)) - diff
    t1,t2 = (lstsmod[0]-0.5*step)*12/n.pi, (lstsmod[-1]+0.5*step)*12/n.pi
else: 
    t1 = (lsts[0]-.5*step)*12/n.pi
    t2 = (lsts[0]+.5*step)*12/n.pi

#lsts = lstsmod

#plot the visibility and phase
p.figure(1)
p.subplot(311); p.plot(lsts,vis.real,label='real vis')
p.subplot(311); p.plot(lsts,vis.imag,label='imag vis')
p.subplot(311); p.plot(lsts,n.abs(vis),label='abs vis')
p.legend()
p.subplot(312); p.plot(lsts, wgts,label='weights')
p.legend()
p.subplot(313); p.plot(lsts,n.angle(vis))


#ifft along the time axis to show the sources.
#First window the data
tw = a.dsp.gen_window(vis.shape[0], window='blackman-harris'); tw.shape += (1,)
tw = a.dsp.gen_window(vis.shape[0], window='blackman-harris'); tw.shape += (1,)
vis_tft = n.fft.fftshift(n.fft.ifft(vis*tw,axis=0),axes=0)
#get fr bins for plotting.
frates = n.fft.fftshift(n.fft.fftfreq(len(vis), 42.8) * 1e3) #mHz
step = frates[1] - frates[0]
f1, f2 = frates[0] - .5*step, frates[-1]+.5*step

#plot fringe rate transform.
#p.figure(3)
#p.subplot(231); p.plot(frates, vis_tft, label='fringe rate transform of vis')
#p.subplot(231); p.plot(frates, n.abs(vis_tft), label='abs')
#:p.title('Stokes I in Fringe rate domain (mHz)')
#p.show()


#Now filter with aggressive fringe rate filter. 
b1 = 1; b2 = 4 #sep24 is the 0,1 spacing.
beam_w_fr = frf_conv.get_beam_w_fr(aa, (b1,b2)) #returns the beam weighted fringe rates in fr domain.

#using a filter that is 3.5785bar hours long = 301 samples.
tbins, firs, frbins, frspace = frf_conv.get_fringe_rate_kernels(beam_w_fr, 42.8, 401)
fr_bins = n.fft.fftshift(n.fft.fftfreq(401, 42.8))

bwfr10 = beam_w_fr.take(chans, axis=0)[0](fr_bins)

firs = firs.take(chans, axis=0).T
frspace= frspace.take(chans, axis=0).T
firs = firs.conj()
#show filters in both spaces.
print firs.shape
print frspace.shape

fq = n.average(freqs)
print fq

#Aarons methof with the sigma fixed
#t = n.arange(-200,200,1) * 42.8
#w = a.dsp.gen_window(t.size, 'none')
#sig = .000288*(fq/.1788)#was found in fr
#cen = .001059 * (fq/.1788)#was found in fr
#firg = n.exp(-(t**2)/(2*(1./(2*n.pi*sig)**2))).astype(n.complex) * w #note that this is in time.
#firg /= firg.sum()
#firg *= n.exp(2j*n.pi*cen*t) # need to flip the sign for backward baselines
#firg = firg.conj()


#firs = firg

vis_frf = n.zeros_like(vis)
wgts_frf = n.zeros_like(vis)
filtered = n.zeros_like(vis)
for ch in xrange(len(chans)):
    #vis_frf[:,ch] = n.convolve(wgts[:,ch]*vis[:,ch], n.conj(firs[ch,:]), mode='same')
    #wgts_frf[:,ch] = n.convolve(wgts[:,ch], n.abs(n.conj(firs[ch,:])), mode='same')
    #filtered[:,ch] = n.where(wgts_frf[:,ch]>0, vis_frf[:,ch]/wgts_frf[:,ch], 1)
    vis_frf[:,ch] = n.convolve(wgts[:,ch]*vis[:,ch], firs[:,ch], mode='same')
    wgts_frf[:,ch] = n.convolve(wgts[:,ch], n.abs(firs[:,ch]), mode='same')
    filtered[:,ch] = n.where(wgts[:,ch]>0, vis_frf[:,ch]/wgts_frf[:,ch], 1)
ww = n.where(wgts_frf>.1, 1, 0)
filtered = ww*filtered


filtered_in_fr = n.fft.fftshift(n.fft.ifft(filtered*tw,axis=0)) #tw is a windowing function.

frbins = frbins*1000 #put in mhz

lsts = lsts*12/n.pi

p.figure(2)
p.suptitle('Convolution method')
p.subplot(241)
p.plot(lsts, vis.real, label='real')
p.plot(lsts, vis.imag, label='imag')
p.plot(lsts, n.abs(vis), label='abs')
p.title('Vis')
p.subplot(242)
p.plot(frates, vis_tft.real, label='real')
p.plot(frates, vis_tft.imag, label='imag')
p.plot(frates, n.abs(vis_tft), label='abs')
p.title('Vis in FR')
p.subplot(243)
p.plot(frbins, frspace.real, label='real')
p.plot(frbins, frspace.imag, label='imag')
p.plot(frbins, n.abs(frspace), label='abs')
p.title('fr filter')
p.subplot(244)
p.plot(frates, filtered_in_fr.real, label='real')
p.plot(frates, filtered_in_fr.imag, label='imag')
p.plot(frates, n.abs(filtered_in_fr), label='abs')
p.title('filtered in fr')


p.subplot(245)
p.plot(frbins, bwfr10.real, label='real')
p.plot(frbins, bwfr10.imag, label='real')
p.plot(frbins, n.abs(bwfr10), label='real')
#p.plot(lsts, vis.real, label='real')
#p.plot(lsts, vis.imag, label='imag')
#p.plot(lsts, n.abs(vis), label='abs')
p.title('Vis')
p.subplot(246)
p.plot(lsts, vis.real, label='real')
p.plot(lsts, vis.imag, label='imag')
p.plot(lsts, n.abs(vis), label='abs')
p.title('Vis')
p.subplot(247)
p.plot(tbins, firs.real, label='real')
p.plot(tbins, firs.imag, label='imag')
p.plot(tbins, n.abs(firs), label='abs')
p.title('time filter')
p.subplot(248)
p.plot(lsts, filtered.real, label='real')
p.plot(lsts, filtered.imag, label='imag')
p.plot(lsts, n.abs(filtered), label='abs')
p.title('filtered in time')
p.savefig('fr_preserved_signal.png', format='png')


p.figure(5)
p.subplot(211);p.plot(frbins, n.abs(frspace))
p.xlabel('Fringe Rate (mHz)')
p.ylabel('amplitude')
p.subplot(212)
p.plot(tbins, n.abs(firs), label='abs')
p.xlabel('Time (sec)')
p.ylabel('amplitude')
p.plot(tbins, firs.real, label='real')
p.plot(tbins, firs.imag, label='imag')
p.suptitle('Fringe Rate filter at 159 MHz')
p.savefig('fr_filter_slice.png', format='png')

#FFT Method
#fft_filter = n.fft.fft(wgts[:,ch]*vis[:,ch]) * firs
#from scipy.signal import fftconvolve
#fft_filter = fftconvolve(wgts[:,0]*vis[:,0], n.squeeze(frspace[:,0]), mode='same')
#fft_filter.shape += (1,)

#import IPython
#IPython.embed()


#p.figure(3)
#p.suptitle('FFT Method')
#p.subplot(241)
#p.plot(lsts, vis.real, label='real')
#p.plot(lsts, vis.imag, label='imag')
#p.plot(lsts, n.abs(vis), label='abs')
#p.title('Vis')
#p.subplot(242)
#p.plot(frates, vis_tft.real, label='real')
#p.plot(frates, vis_tft.imag, label='imag')
#p.plot(frates, n.abs(vis_tft), label='abs')
#p.title('Vis in FR')
#p.subplot(243)
#p.plot(frbins, frspace.real, label='real')
#p.plot(frbins, frspace.imag, label='imag')
#p.plot(frbins, n.abs(frspace), label='abs')
#p.title('fr filter')
#p.subplot(244)
#p.plot(frates, n.fft.ifft(fft_filter*tw).real, label='real')
#p.plot(frates, n.fft.ifft(fft_filter*tw).imag, label='imag')
#p.plot(frates, n.abs(n.fft.ifft(fft_filter*tw)), label='abs')
#p.title('filtered in fr')
#
#
#p.subplot(245)
#p.plot(lsts, vis.real, label='real')
#p.plot(lsts, vis.imag, label='imag')
#p.plot(lsts, n.abs(vis), label='abs')
#p.title('Vis')
#p.subplot(246)
#p.plot(lsts, vis.real, label='real')
#p.plot(lsts, vis.imag, label='imag')
#p.plot(lsts, n.abs(vis), label='abs')
#p.title('Vis')
#p.subplot(247)
#p.plot(tbins, firs.real, label='real')
#p.plot(tbins, firs.imag, label='imag')
#p.plot(tbins, n.abs(firs), label='abs')
#p.title('time filter')
#p.subplot(248)
#p.plot(lsts, filtered.real, label='real')
#p.plot(lsts, filtered.imag, label='imag')
#p.plot(lsts, n.abs(filtered), label='abs')
#p.title('filtered in time')


#do fft filter.
#t, firs, frbins, frspace = frf_conv.get_fringe_rate_kernels(beam_w_fr, 42.8, len(vis))
#frspace = frspace[chans][0]
#p.figure(7);p.plot(frates, n.fft.fft(wgts[:,ch]*vis[:,ch] * frspace))
#p.figure(7);p.plot(frates, n.abs(n.fft.fft(wgts[:,ch]*vis[:,ch] * frspace)))
#p.figure(8);p.plot(frbins, frspace)
#p.figure(8);p.plot(frbins, n.abs(frspace))
p.figure(4,figsize=(14,6))
p.suptitle('Filter Steps')
p.subplot(131)
p.plot(frates, vis_tft.real, label='real')
p.plot(frates, vis_tft.imag, label='imag')
p.plot(frates, n.abs(vis_tft), label='abs')
p.xlim(-.6,1.5)
p.xlabel('Fringe Rate  (mHz)')
p.ylabel('amplitude')
p.title('Vis in FR')
p.subplot(132)
p.plot(frbins, frspace.real, label='real')
p.plot(frbins, frspace.imag, label='imag')
p.plot(frbins, n.abs(frspace), label='abs')
p.xlim(-.6,1.5)
p.xlabel('LST  (hrs)')
p.ylabel('amplitude')
p.title('fr filter')
p.xlabel('Fringe Rate  (mHz)')
p.ylabel('amplitude')
p.subplot(133)
p.plot(frates, filtered_in_fr.real, label='real')
p.plot(frates, filtered_in_fr.imag, label='imag')
p.plot(frates, n.abs(filtered_in_fr), label='abs')
p.xlim(-.6,1.5)
p.xlabel('Fringe Rate  (mHz)')
p.ylabel('amplitude')
p.title('filtered in fr')


#p.subplot(234)
#p.plot(lsts, vis.real, ',', label='real')
#p.plot(lsts, vis.imag, ',', label='imag')
#p.plot(lsts, n.abs(vis), ',', label='abs')
#p.xlabel('LST  (hrs)')
#p.ylabel('amplitude')
#p.title('Vis')
#p.subplot(235)
#p.plot(tbins, firs.real, label='real')
#p.plot(tbins, firs.imag, label='imag')
#p.plot(tbins, n.abs(firs), label='abs')
#p.xlabel('time  (sec)')
#p.ylabel('amplitude')
#p.title('Convolution Kernel')
#p.subplot(236)
#p.plot(lsts, filtered.real, ',', label='real')
#p.plot(lsts, filtered.imag, ',', label='imag')
#p.plot(lsts, n.abs(filtered), ',', label='abs')
#p.xlabel('LST  (hrs)')
#p.ylabel('amplitude')
#p.title('filtered in time')

p.savefig('fr_preserved_signal.png', format='png')

#p.show()






