#feed this script uv files to make noise from
#Does not require sample number since deleteing simuations after using them.

for DAY in even odd; do 
    cd signal/sim/${DAY}
    echo rm -rf lst*.uvALG_signal
    rm -rf lst*.uvALG_signal
    echo rm -rf lst*.uvALG_signalL
    rm -rf lst*.uvALG_signalL
    echo make_noise_uv.py --signal ../../../data/${DAY}/lst*.uvALG
    make_noise_uv.py --signal ../../../data/${DAY}/lst*.uvALG
    echo fringe_rate_filter_convolution.py -C psa6240_v003 --sep=0,1 *.uvALG_signal
    fringe_rate_filter_convolution.py -C psa6240_v003 --sep=0,1 *.uvALG_signal
#    echo cd ../../..
    cd ../../..
done

