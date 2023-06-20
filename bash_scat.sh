LANG=en_US
npini=6
npadd=20
p0=0.000005
prep=50
np=$(echo "$npini + $npadd" | bc -l)
del=0.1
al=0.1
sig=0.005
tmax=1500
nmode=1200
wc=1
wmax=3
x0=-700
xmin=400
nval=0.5
tref=3
me=0.0000001
kmin=0.06025
kmax=0.20
kint=0.001
knumb=$(echo "($kmax-$kmin)/$kint" | bc -l)
x0char=$(python -c "print -1*$x0/1000")
echo $knumb

dt=0.05
dtchar=$(python -c "print $dt*100")


MPOL_DIR_SCRATCH=/scratch/gheeraert/mpoldyn_cluster_27
paramchar=scat_n$nval\x$x0char\nm$nmode\np$npini\_$npadd\_tref$tref\_al$al\sig$sig\dt$dt\t$tmax\p$p0
paramchar_assembled=scat_n$nval\x$x0char\nm$nmode\np$npini\_$npadd\_al$al\sig$sig\dt$dt\t$tmax
WORK_DIR=$MPOL_DIR_SCRATCH/data/$paramchar
mkdir $WORK_DIR
mkdir $WORK_DIR/data
log=$WORK_DIR/data/log_$paramchar

SCAT_DIR=$MPOL_DIR_SCRATCH/data/scat_folder
mkdir $SCAT_DIR

cp -r $MPOL_DIR_SCRATCH/gs_data $WORK_DIR

for i in $(seq 0 $knumb);
do
echo $i

kvalue=$(echo "$kmin + $kint*$i" | bc -l)
echo "kvalue="
echo $kvalue
kvalformated=$(printf '%6.4f\n' $kvalue)
echo $(printf '%6.4f\n' $kvalue)
echo "kvalueformated="
echo $kvalformated
  log=$WORK_DIR/data/log_$paramchar\_$kvalformated
  jobfilename=$WORK_DIR/RT_np$np\_n$nval\_$kvalue

  cp $MPOL_DIR_SCRATCH/jobfile $jobfilename
  echo "cd $WORK_DIR" >> $jobfilename

  echo "$MPOL_DIR_SCRATCH/mpol -n $nval -x0 $x0 -me $me -sigma $sig -prep $prep -tmax $tmax -k0 $kvalue -xmin $xmin -npini $npini -npadd  $npadd -p0 $p0 -tref $tref -del $del -nm $nmode -wc $wc -wmax $wmax -al $al -dt $dt > $log " >> $jobfilename
  echo "$MPOL_DIR_SCRATCH/mpol -n $nval -x0 $x0 -me $me -sigma $sig -prep $prep -tmax $tmax -k0 $kvalue -xmin $xmin -npini $npini -npadd  $npadd -p0 $p0 -tref $tref -del $del -nm $nmode -wc $wc -wmax $wmax -al $al -dt $dt > $log "

  qsub $jobfilename

  echo  " awk 'NR==1 {print  "'$1'", "'$2'", "'$3'", "'$4'", "'$5'","'$6'","'$7'","'$8'","'$9'"}' ../$paramchar/file_$kvalformated.d >> ./$paramchar_assembled.d " >>  $SCAT_DIR/assemble.sh

done

paramchar_sorted=scat_n$nval\x$x0char\nm$nmode\np$npini\_$npadd\_al$al\sig$sig\dt$dt\t$tmax\p$prep\_sorted
echo  "sort -k1 -n $paramchar_assembled.d > $paramchar_sorted.d" >>  $SCAT_DIR/assemble.sh
echo  "rm $paramchar_assembled.d" >>  $SCAT_DIR/assemble.sh

#gnumb=$(echo "($gmax-$gmin)/$gint" | bc -l)
#gnumb=$(python -c "print int($gnumb)")
#
#
#for i in $(seq 0 $gnumb);
#do
#gvalue=$(echo "$gmin + $gint*$i" | bc -l)
#
#cp $WDIR/jobfile $WDIR2/jf_g$gvalue\_nS$nSteps\_v$vCoul\_i$i
#echo "$WDIR/exdy -L $nSites -h $hval -m $mval -w $wval -v $vCoul -dt 0.01 -nS $nSteps -m $mval -g $gvalue > $WDIR2/tmp_v$vCoul\_nS$nSteps\_L$nSites\_m$mval\_w$wval\_g$gvalue\.d" >> $WDIR2/jf_g$gvalue\_nS$nSteps\_v$vCoul\_i$i
#
#qsub jf_g$gvalue\_nS$nSteps\_v$vCoul\_i$i 
#
#dond:w

