# ==============================================================================
#  FILE: bash_1scat.sh
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
npini=6
# prep: Initialisation mode passed to the executable (-prep). In this codebase, 50/51 select
# coherent or single-photon pulses on top of the dressed ground state.
prep=50
# del: Bare qubit splitting Δ (passed as -del).
del=0.1
# al: Dimensionless coupling α (passed as -al).
al=0.1
# tmax: Final simulation time.
tmax=2600
# nmode: Number of discretised bosonic modes used to approximate the continuum.
nmode=2500
# wc: High-frequency cutoff ωc (sets the exponential cutoff of the spectral density).
wc=1
# wmax: Maximum frequency retained in the discretised band.
wmax=3
# x0: Initial centre position of the incident wavepacket (typically negative, to start left of the
# scatterer).
x0=-1400
# xmin: Spatial cutoff used when computing Fourier transforms / field profiles (filters far-field
# region).
xmin=600
# nval: Average photon number / pulse intensity parameter used for coherent-state input (-n).
nval=0.5
# merr: Error threshold used to trigger basis enlargement (-me).
merr=0.0000001
x0char=$(python -c "print -1*$x0/1000")
# mde: Maximum allowed energy drift; used as an integration/basis-stability guard (-mde).
mde=0.00000001
# p0: Initial amplitude assigned to newly-added coherent components (-p0).
p0=0.000001

#k0=0.2075
# k0: Central wavenumber/frequency of the incident pulse (near resonance).
k0=0.166
# sig: Spectral width of the incident pulse (Gaussian in k-space). Passed as -sigma.
sig=0.005
trefmin=0.95
trefmax=3.0
trefint=0.1
trefnumb=$(echo "($trefmax-$trefmin)/$trefint" | bc -l)

# dt: Base integration time step (-dt).
dt=0.04
dtchar=$(python -c "print $dt*100")

# MPOL_DIR_SCRATCH: Working directory on the cluster containing the executable and templates.
MPOL_DIR_SCRATCH=/scratch/gheeraert/mpoldyn_cluster_116

for i in $(seq 8 8);
do

  npadd=$(echo "$i*3" | bc -l)
  np=$(echo "$npini+$npadd" | bc -l)
  echo i=$i, np=$np

  for j in $(seq 0 $trefnumb);
  do
	 tref=$(echo "$trefmin + $trefint*$j" | bc -l)
	 echo  ---j=$j, tref=$tref

	 jobfilename=$MPOL_DIR_SCRATCH/$np\_$k0\_$tref\_$al
	 cp $MPOL_DIR_SCRATCH/jobfile $jobfilename
	 paramchar=np$np\_nm$nmode\_al$al\_del$del\_k$k0\_x$x0\_n$nval\_sig$sig\_t$tmax\_tr$tref\_me$merr\_dt$dtchar
	 log=$MPOL_DIR_SCRATCH/data/log_$paramchar.d
	 echo "cd $MPOL_DIR_SCRATCH" >> $jobfilename
	 echo "module load intel/14.0.1" >> $jobfilename
	 echo "$MPOL_DIR_SCRATCH/mpol -n $nval -x0 $x0 -me $merr -sigma $sig -prep $prep -tmax $tmax -k0 $k0 -xmin $xmin -npini $npini -npadd  $npadd -tref $tref -del $del -nm $nmode -wc $wc -wmax $wmax -al $al -dt $dt -mde $mde -p0 $p0 > $log " >> $jobfilename
	 echo "$MPOL_DIR_SCRATCH/mpol -n $nval -x0 $x0 -me $merr -sigma $sig -prep $prep -tmax $tmax -k0 $k0 -xmin $xmin -npini $npini -npadd  $npadd -tref $tref -del $del -nm $nmode -wc $wc -wmax $wmax -al $al -dt $dt -mde $mde -p0 $p0 > $log "
	 qsub $jobfilename
	 rm $jobfilename

  done

done
