#! /usr/bin/env python
import numpy as n
import capo as C 
import optparse,sys,glob

import pylab as p


o = optparse.OptionParser()
o.add_option('--plot', action='store_true', 
    help='plot stuff')
opts,args = o.parse_args(sys.argv[1:])
PLOT=opts.plot


#load reference data for signal loss, both xc^-1Qc^-1x, and xQx.
data_file_pk = n.load('../plots/data/signalloss/level_0000/pspec.npz')
data_file_nocov = n.load('../plots/data/signalloss/level_0000/nocov_pspec.npz')

pk_ref = data_file_pk['pk']
pk_ref_err = data_file_pk['err']

nocov_ref = data_file_nocov['pk']
nocov_ref_err = data_file_nocov['err']

#k_parallel's
kpl = data_file_pk['kpl']

#want to eventually monte carlo when we get enough simulations.
# Hence we want the data products: 
#    1) P_c(x+xs) - P_c(x)
#    2) P_i(x+xs) - P_i(x) 
#note that P_c(x) and P_i(x) do not change.Want 1/2.

#dir structure: v000{i}. - each version is a different signal level input. 
# each version has multiple instanstiations of the simulatied signal.

#pks, nocovs, ratios structure:
#   ratios[level] => level is inpput signal level. For each level, there is a a list of the samples => shape = (nsamples, nkpls) 
pks = {}
#pk_errs = {}
nocovs = {}
#nocov_errs = {}
ratios = {}
for level in args:
    if level=='level_0000': continue
    print 'Reading files from ', level
    if not level in pks: 
        pks[level] = []
        nocovs[level] = []
        ratios[level] = []
        #pk_errs[level] = {}
        #nocov_errs[level] = {}
    #for sample in glob.glob(level+'/sample_[0-9]')+glob.glob(level+'/sample_[1-7][0-9]'):
    for sample in glob.glob(level+'/sample_1[0-3][0-9]'):
        pkfile = n.load(sample+'/pspec.npz')
        nocovfile = n.load(sample+'/nocov_pspec.npz')

        pks[level].append(pkfile['pk'])
        nocovs[level].append(nocovfile['pk'])
    
        #pk_errs[level][sample] = pkfile['err']
        #nocov_errs[level][sample] = nocovfile['err']

#        ratios[level].append((pkfile['pk'] - pk_ref)/(nocovfile['pk'] - nocov_ref))



avg_ratios = []
var_ratios = []
avg_pks = []
avg_nocov = []
var_pks = []
var_nocov = []
lkeys = []
#p.plot(n.transpose(ratios[ratios.keys()[0]]))
#p.show()
for level in pks.keys():
    pks[level] = n.array(pks[level])
    nocovs[level] = n.array(nocovs[level])
#    ratios[level] = n.array(ratios[level])
#    print ratios[level].shape
    print pks[level].shape
#    avg_ratios.append(n.mean(ratios[level], axis=0))
#    var_ratios.append(n.var(ratios[level], axis=0))
    avg_pks.append(n.mean(pks[level], axis=0))
    var_pks.append(n.var(pks[level], axis=0))
    avg_nocov.append(n.mean(nocovs[level], axis=0))
    var_nocov.append(n.var(nocovs[level], axis=0))
    lkeys.append(level.split('_')[-1])
if PLOT:
    for pk,nocov in zip(avg_pks,avg_nocov):
        p.figure(1)
        p.plot(avg_pks[-1], label=level)
        p.figure(2)
        p.plot(avg_nocov[-1], label=level)
    p.figure(1)
    p.plot(pk_ref, '--')
    p.legend()
    p.figure(2)
    p.plot(nocov_ref, '--')
    p.legend()
    p.show()

avg_pks = n.array(avg_pks)
avg_nocov= n.array(avg_nocov)
var_pks = n.array(var_pks)
var_nocov = n.array(var_nocov)
#avg_ratios = n.array(avg_ratios)
#var_ratios = n.array(var_ratios)

print 'Plotting'

#for k in range(9):
    #p.errorbar(n.arange(avg_ratios.shape[1]), avg_ratios[k], yerr = var_ratios[k])
#    p.plot(n.arange(avg_ratios.shape[1]), avg_ratios[k], 'o')

#p.show()


loss_ratio = []
#get_k = []
for k, level in enumerate(pks.keys()):
    print k, level
    loss_ratio.append( (avg_nocov[k] - nocov_ref)/(avg_pks[k] - pk_ref) )
#    p.errorbar(loss_ratio[-1][15], label=level)
#    get_k.append(loss_ratio[-1][15])
loss_ratio = n.array(loss_ratio)

if PLOT:
    for k, level in enumerate(pks.keys()):
        p.figure(1)
        p.plot( (avg_pks[k] - pk_ref), label=level)
        p.figure(2)
        p.plot( (avg_nocov[k] - nocov_ref), label=level)
        p.figure(3)
        p.plot( loss_ratio[k], label=level)
    p.figure(1)
    p.legend()
    p.figure(2)
    p.legend()
    p.figure(3)
    p.legend()
    p.show()


#for i,d in enumerate(avg_ratios):
lkeys = n.array([float(k) for k in lkeys])
argkeys = n.argsort(lkeys)
print argkeys
loss_ratio = loss_ratio[argkeys]
print loss_ratio.shape
#p.plot(n.arange(avg_ratios.shape[0]), n.mean(avg_ratios,axis=1))#, label=lkeys[i])
signal_levels = lkeys[argkeys]
average_these = n.ones(21)
average_these[10-4:10] = 0
average_these[10:10+5] = 0
p.plot(n.log(signal_levels), n.average(loss_ratio,axis=1,weights=average_these))
#p.errorbar(n.log10(signal_levels), n.average(loss_ratio,axis=1,weights=average_these), yerr = n.var(loss_ratio, axis=1))#, label=lkeys[i])
p.legend()
p.show()
