# ==============================================================================
#  FILE: bash_1scat.sh
#
#  PURPOSE
#    Submission / sweep script for running the Fortran executable (mpol) on a
#    cluster scheduler (qsub). The generated data files are then used for plotting
#    scattering observables (reflection / transmission) and photon statistics.
#    
# ==============================================================================
#
npini=6         # npini: Initial number of coherent-state components ("polarons") in the variational ansatz.
prep=50         # prep: Initialisation mode passed to the executable (-prep).
del=0.1         # del: Bare qubit splitting Δ (passed as -del).
al=0.1          # al: Dimensionless coupling α (passed as -al).
tmax=2600       # tmax: Final simulation time.
nmode=2500      # nmode: Number of discretised bosonic modes used to approximate the continuum.
wc=1            # wc: High-frequency cutoff ωc (sets the exponential cutoff of the spectral density).
wmax=3          # wmax: Maximum frequency retained in the discretised band.
x0=-1400        # x0: Initial centre position of the incident wavepacket
xmin=600        # xmin: Spatial cutoff used when computing Fourier transforms / field profiles
nval=0.5        # nval: Average photon number / pulse intensity parameter used for coherent-state input (-n).
merr=0.0000001  # merr: Error threshold used to trigger basis enlargement (-me).
x0char=$(python -c "print -1*$x0/1000")
mde=0.00000001  # mde: Maximum allowed energy drift; used as an integration/basis-stability guard (-mde).
p0=0.000001     # p0: Initial amplitude assigned to newly-added coherent components (-p0).
k0=0.166        # k0: Central wavenumber/frequency of the incident pulse (near resonance).
sig=0.005       # sig: Spectral width of the incident pulse (Gaussian in k-space). Passed as -sigma.
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
