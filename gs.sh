npini=5
del=0.1
al=0.1
nmode=300
wc=1000
wmax=1

WDIR_gs=/home/gheeraert/gs_mpol
WDIR=/home/gheeraert/mpoldyn_cluster

gs_jobfile=$WDIR/gs_data/jf_np$npini
cp $WDIR/jobfile $gs_jobfile

echo "cd $WDIR/gs_data" >> $gs_jobfile
echo "$WDIR_gs/mpol -npol $npini -nk $nmode -wc $wc -wmax $wmax -alpha $al -delta $del " >> $gs_jobfile

qsub $gs_jobfile

rm $gs_jobfile

