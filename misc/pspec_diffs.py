#! /usr/bin/env python
import numpy as n
import capo as C 
import optparse,sys,glob

import pylab as p


o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

#load reference data for signal loss, both xc^-1Qc^-1x, and xQx.
data_file_pk = n.load('level_0000/pspec.npz')
data_file_nocov = n.load('level_0000/nocov_pspec.npz')

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
    for sample in glob.glob(level+'/sample_*'):
        pkfile = n.load(sample+'/pspec.npz')
        nocovfile = n.load(sample+'/nocov_pspec.npz')

        pks[level].append(pkfile['pk'])
        nocovs[level].append(nocovfile['pk'])
    
        #pk_errs[level][sample] = pkfile['err']
        #nocov_errs[level][sample] = nocovfile['err']

        ratios[level].append((pkfile['pk'] - pk_ref)/(nocovfile['pk'] - nocov_ref))



avg_ratios = []
var_ratios = []
lkeys = []
#p.plot(n.transpose(ratios[ratios.keys()[0]]))
#p.show()
for level in ratios.keys():
    ratios[level] = n.array(ratios[level])
    print ratios[level].shape
    avg_ratios.append(n.mean(ratios[level], axis=0))
    var_ratios.append(n.var(ratios[level], axis=0))
    lkeys.append(level.split('_')[-1])

avg_ratios = n.array(avg_ratios)
var_ratios = n.array(var_ratios)

print 'Plotting'
print avg_ratios.shape

#for k in range(9):
    #p.errorbar(n.arange(avg_ratios.shape[1]), avg_ratios[k], yerr = var_ratios[k])
#    p.plot(n.arange(avg_ratios.shape[1]), avg_ratios[k], 'o')

#p.show()



#for i,d in enumerate(avg_ratios):
lkeys = n.array([float(k) for k in lkeys])
argkeys = n.argsort(lkeys)
print argkeys
avg_ratios = avg_ratios[argkeys]
print avg_ratios.shape
#p.plot(n.arange(avg_ratios.shape[0]), n.mean(avg_ratios,axis=1))#, label=lkeys[i])
signal_levels = lkeys[argkeys]
p.errorbar(n.log10(signal_levels), n.mean(avg_ratios,axis=1), yerr = n.var(avg_ratios, axis=1))#, label=lkeys[i])
p.legend()
p.show()
