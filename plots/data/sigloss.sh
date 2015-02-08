#! /usr/bin/env bash

for inject in `python -c "import numpy; print ' '.join(map(str, numpy.logspace(-2,2,40)))"` ; do
    echo $inject
    python ~/pkgs/capo/arp/scripts/pspec_cov_v002_sigloss.py -c 95_115 --window=none -C psa6240_v003 -b 40 -i $inject
    mkdir sigloss/inject_${inject}
    echo "python ~/pkgs/capo/arp/scripts/pspec_cov_v002_sigloss.py -c 95_115 --window=none -C psa6240_v003 -i $inject" > sigloss/inject_${inject}/notes.txt
    mv pspec_boot*.npz sigloss/inject_${inject}
done
