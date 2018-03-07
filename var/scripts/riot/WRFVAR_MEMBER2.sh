#!/bin/bash
SECONDS=0
#========================================================
#WRFVAR_MEMBER.sh is used to run a single ensemble member
#========================================================
iSAMP="$1"; npiens="$2";

ii=$iSAMP
if [ $iSAMP -lt 10 ]; then ii="0"$ii; fi
if [ $iSAMP -lt 100 ]; then ii="0"$ii; fi
if [ $iSAMP -lt 1000 ]; then ii="0"$ii; fi

cd $MEMBERPREFIX$ii

cp $CWD/namelist.input ./
ex -c :"/ensmember" +:":s/=.*/=$iSAMP,/" +:wq namelist.input

# Multiply A^T * A * q_i [iSAMP<=NSAMP] or A^T * R^-1/2 * [h(x + dx) - y + dy] [iSAMP==NSAMP+1]
# A = R^-1/2 * H * L
# dx and dy are only used in Block Lanczos (rand_type==3)
# Note: The implementation through mpirun or mpiexec is 
#       unique for your cluster and MPI implementation

# DEFAULT - WORKS FOR:
#  1a NASA Pleiades for SGI MPT
#  1b NOAA Theia for Intel MPI
#  *Please add more as verified or other solutions developed

export PBS_NODEFILE=$(pwd)/hostlist #.$ii

#BACKG_STRING="< /dev/null &> out.run &"
BACKG_STRING="< /dev/null &> out.run"

mpistring="$MPICALL $DEBUGSTR -np $npiens $RIOT_EXECUTABLE $BACKG_STRING"

# 2 - Potentially generic solution for Intel MPI implmentations
#      (may need earlier calls to mpdallexit and mpdboot)
#mpistring="$MPICALL -np $npiens -machinefile $(pwd)/hostlist $RIOT_EXECUTABLE $BACKG_STRING"

# 3 - Potentially generic solution for Open MPI implmentations
#      (may need earlier calls to mpdallexit and mpdboot)
#mpistring="$MPICALL -np $npiens --hostfile $(pwd)/hostlist $RIOT_EXECUTABLE $BACKG_STRING"

#Use these if $BACKGR_STRING does not place mpirun in background (no trailing "&"):
TEMPTIME=$SECONDS
echo "JOB $iSAMP; $mpistring; $(date); $SECONDS sec."
eval "$mpistring"
TEMPTIME=$((SECONDS-TEMPTIME))

hr_temp0=$(($TEMPTIME / 3600))
if [ $hr_temp0 -lt 10 ]; then hr_temp0="0"$hr_temp0; fi
min_temp0=$((($TEMPTIME / 60) % 60))
if [ $min_temp0 -lt 10 ]; then min_temp0="0"$min_temp0; fi
sec_temp0=$(($TEMPTIME % 60))
if [ $sec_temp0 -lt 10 ]; then sec_temp0="0"$sec_temp0; fi

##   #Use these if $BACKGR_STRING places mpirun in background (trailing "&", no timing possible):
##   eval "$mpistring"; wait_pids+=($!)
##   echo "JOB $iSAMP; wait ${wait_pids[@]}"
##   wait "${wait_pids[@]}"

#NOTE: Redirected input ($BACKG_STRING) necessary to run in background and ensures clean log files

mpireturn=$?

#Reset PBS_NODEFILE
export PBS_NODEFILE=$PBSNODE0

hr0=$(($SECONDS / 3600))
if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
min0=$((($SECONDS / 60) % 60))
if [ $min0 -lt 10 ]; then min0="0"$min0; fi
sec0=$(($SECONDS % 60))
if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
echo "WRFVAR_MEMBER.$ii time, TOTAL-$hr0:$min0:$sec0; "$'\n'"WRFDA-$hr_temp0:$min_temp0:$sec_temp0; "$'\n'"w/ WRFDA return value: $mpireturn"

exit "$mpireturn"

