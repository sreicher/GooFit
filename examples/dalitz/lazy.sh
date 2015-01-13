rm dalitz
#gmake clean && gmake
gmake clean && gmake CUDALOCATION=/afs/hep.man.ac.uk/u/reichert/OMP_GooFit/GooFit/fakecuda/ TARGET_OMP=1
./dalitz Toy_noPhotos_isobar_x0.005_y0.005.txt SL_KsPiPi_BGun_isobar_x0.005_y0.005_mptp_no_time-dep_ 
#./dalitz SL_KsPiPi_BGun_isobar_x0.005_y0.005_mptp_no_time-dep_788.txt
