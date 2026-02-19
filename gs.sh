# ==============================================================================
#  FILE: gs.sh
#
#  PURPOSE
#    Submission / sweep script for running the Fortran executable (mpol) on a
#    cluster scheduler (qsub). The generated data files are then used for plotting
#    scattering observables (reflection / transmission) and photon statistics.
#
#  CONTEXT (paper)
#    The parameters correspond to the spin-boson / waveguide-QED model and the
#    drive/wavepacket protocol described in Gheeraert et al., PRA 98, 043816 (2018).
#
#  Generated on: 2026-02-19
# ==============================================================================
#
# npini: Initial number of coherent-state components ("polarons") in the variational ansatz.
npini=5
# del: Bare qubit splitting Δ (passed as -del).
del=0.1
# al: Dimensionless coupling α (passed as -al).
al=0.1
# nmode: Number of discretised bosonic modes used to approximate the continuum.
nmode=300
# wc: High-frequency cutoff ωc (sets the exponential cutoff of the spectral density).
wc=1000
# wmax: Maximum frequency retained in the discretised band.
wmax=1

WDIR_gs=/home/gheeraert/gs_mpol
WDIR=/home/gheeraert/mpoldyn_cluster

gs_jobfile=$WDIR/gs_data/jf_np$npini
cp $WDIR/jobfile $gs_jobfile

echo "cd $WDIR/gs_data" >> $gs_jobfile
echo "$WDIR_gs/mpol -npol $npini -nk $nmode -wc $wc -wmax $wmax -alpha $al -delta $del " >> $gs_jobfile

qsub $gs_jobfile

rm $gs_jobfile

