import frf_conv
import aipy as a, numpy as n, pylab as p
import capo as C

#get antenna array for frf filters.
aa = a.cal.get_aa('psa6240_v003', .1/203, .1, 203)
#channel range
chans = n.arange(90,130)
#chans = n.arange(110,150)
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
vis = n.concatenate([vis[-200:], vis[:1460]], axis=0) #want in time order for fr filter
wgts = n.concatenate([wgts[-200:], wgts[:1460]], axis=0) #want in time order for frf
lsts = n.concatenate([lsts[-200:], lsts[:1460]], axis=0) #want in time order for frf

lst_rng = (22,5)
lst_mask = n.logical_and(lsts>=3, lsts<=5)

vis = vis[lst_mask]
wgts = wgts[lst_mask]
lsts = lsts[lst_mask]
print vis.shape



#for plotting the correct lsts.
step = lsts[1] - lsts[0]
if lsts[0] > lsts[-1]:#wrapped
    diff = 2*n.pi - lsts[0]
    lstsmod = ((lsts + diff)%(2*n.pi)) - diff
    t1,t2 = (lstsmod[0]-0.5*step)*12/n.pi, (lstsmod[-1]+0.5*step)*12/n.pi
else: 
    t1 = (lsts[0]-.5*step)*12/n.pi
    t2 = (lsts[0]+.5*step)*12/n.pi

#plot the visibility and phase
p.figure(1);C.arp.waterfall(vis, extent=(chans[0],chans[-1],t2,t1)); p.colorbar(shrink=.5)
p.figure(2);C.arp.waterfall(vis, mode='phs', extent=(chans[0],chans[-1],t2,t1)); p.colorbar(shrink=.5)

#ifft along the time axis to show the sources.
#First window the data
tw = a.dsp.gen_window(vis.shape[0], window='blackman-harris'); tw.shape += (1,)
vis_tft = n.fft.fftshift(n.fft.ifft(vis*tw,axis=0),axes=0)
vis_del = n.fft.fftshift(n.fft.fft(vis,axis=1),axes=1)
#get fr bins for plotting.
frates = n.fft.fftshift(n.fft.fftfreq(len(vis), 42.8) * 1e3) #mHz
step = frates[1] - frates[0]
f1, f2 = frates[0] - .5*step, frates[-1]+.5*step

#plot fringe rate transform.
p.figure(3);C.arp.waterfall(vis_tft, extent=(chans[0],chans[-1],f2,f1), drng=5);p.colorbar(shrink=.5)
p.title('Stokes I in Fringe rate domain (mHz)')
p.figure(4);C.arp.waterfall(vis_del);p.colorbar(shrink=.5)
#p.show()


#Now filter with aggressive fringe rate filter. 
b1 = 1; b2 = 4 #sep24 is the 0,1 spacing.
beam_w_fr = frf_conv.get_beam_w_fr(aa, (b1,b2)) #returns the beam weighted fringe rates in fr domain.

test = beam_w_fr[110]
#using a filter that is 3.5785bar hours long = 301 samples.
t, firs, frbins, frspace = frf_conv.get_fringe_rate_kernels(beam_w_fr, 42.8, 51)
fr_bins = n.fft.fftshift(n.fft.fftfreq(51, 42.8))

#show filters in both spaces.
print firs.shape
print frspace.shape
p.figure(20)
C.arp.waterfall(firs,extent=(t[0],t[-1],0,203),mx=0,drng=5);p.colorbar(shrink=.5)
p.figure(21)
C.arp.waterfall(frspace,extent=(f1,f2,203,0), mx=0,drng=4.5);p.colorbar(shrink=.5)


firs = firs.take(chans, axis=0)
fq = n.average(freqs)
print fq

#firs = firs[-1]
#firs = firs.conj()
#p.figure(10); #p.plot(fr_bins, test(fr_bins))
##p.plot(fr_bins, n.abs(frspace[110]),'*')
#p.plot(1000*fr_bins, n.abs(n.fft.fft(firs)),'*')
#p.figure(5);
#p.plot(firs.real)
#p.plot(firs.imag)
#p.plot(n.abs(firs))
#firs=firs.conj()


vis_frf = n.zeros_like(vis)
wgts_frf = n.zeros_like(vis)
filtered = n.zeros_like(vis)
for ch in xrange(len(chans)):
    vis_frf[:,ch] = n.convolve(wgts[:,ch]*vis[:,ch], n.conj(firs[ch,:]), mode='same')
    wgts_frf[:,ch] = n.convolve(wgts[:,ch], n.abs(n.conj(firs[ch,:])), mode='same')
    filtered[:,ch] = n.where(wgts_frf[:,ch]>0, vis_frf[:,ch]/wgts_frf[:,ch], 1)
    #vis_frf[:,ch] = n.convolve(wgts[:,ch]*vis[:,ch], firs, mode='same')
    #wgts_frf[:,ch] = n.convolve(wgts[:,ch], n.abs(firs), mode='same')
    #filtered[:,ch] = n.where(wgts[:,ch]>0, vis_frf[:,ch]/wgts_frf[:,ch], 1)
ww = n.where(wgts_frf>.1, 1, 0)
filtered = ww*filtered




#p.figure(4);C.arp.waterfall(n.fft.fftshift(n.fft.ifft(filtered,axis=0), axes=0)); p.colorbar(shrink=.5)
#import IPython
#IPython.embed()

p.figure(6);C.arp.waterfall(n.fft.fftshift(n.fft.ifft(filtered*tw,axis=0),axes=0), extent=(chans[0],chans[-1],f2,f1), drng=5); p.colorbar(shrink=.5)
p.figure(7);C.arp.waterfall(n.fft.fftshift(n.fft.fft(filtered,axis=1),axes=1)); p.colorbar(shrink=.5)
p.figure(8);C.arp.waterfall(filtered,mode='phs'); p.colorbar(shrink=.5)
p.figure(9);C.arp.waterfall(filtered); p.colorbar(shrink=.5)
p.show()






