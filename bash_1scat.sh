npini=6
prep=50
del=0.1
al=0.1
tmax=2600
nmode=2500
wc=1
wmax=3
x0=-1400
xmin=600
nval=0.5
merr=0.0000001
x0char=$(python -c "print -1*$x0/1000")
mde=0.00000001
p0=0.000001

#k0=0.2075
k0=0.166
sig=0.005
trefmin=0.95
trefmax=3.0
trefint=0.1
trefnumb=$(echo "($trefmax-$trefmin)/$trefint" | bc -l)

dt=0.04
dtchar=$(python -c "print $dt*100")

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
