#!/bin/bash
SECONDS=0
#========================================================
#WRFVAR_MEMBER.sh is used to run a single ensemble member
#========================================================
iSAMP="$1"; npiens="$2"; it="$3"; rand_stage="$4"; innerit="$5"

## Alternative compatible with GNU parallel
#   NSAMP="$2"
#   npiens=$nproc_local
#   if [ $1 -eq $((NSAMP+1)) ]; then
#      npiens=$nproc_local_grad
#   fi

ii=$iSAMP
if [ $iSAMP -lt 10 ]; then ii="0"$ii; fi
if [ $iSAMP -lt 100 ]; then ii="0"$ii; fi
if [ $iSAMP -lt 1000 ]; then ii="0"$ii; fi

cd ../run.$ii

it0=$it
if [ $it -lt 10 ]; then it0="0"$it0; fi

it0_last=$((it-1))
if [ $it0_last -lt 10 ]; then it0_last="0"$it0_last; fi

innerit0=$innerit
if [ $innerit -lt 10 ]; then innerit0="0"$innerit0; fi
if [ $innerit -lt 100 ]; then innerit0="0"$innerit0; fi
if [ $innerit -lt 1000 ]; then innerit0="0"$innerit0; fi

innerit0_last=$((innerit-1))
if [ $((innerit-1)) -lt 10 ]; then innerit0_last="0"$innerit0_last; fi
if [ $((innerit-1)) -lt 100 ]; then innerit0_last="0"$innerit0_last; fi
if [ $((innerit-1)) -lt 1000 ]; then innerit0_last="0"$innerit0_last; fi

# Archive old log files
if [ $innerit -eq 1 ]; then
   if [ $rand_stage -eq 1 ]; then
      if [ $it -gt 1 ] && \
         ([ $RIOT_RESTART -eq 0 ] || [ $RIOT_RESTART -eq 2 ]); then
         mkdir oldrsl_$it0_last
         mv rsl.* oldrsl_$it0_last
      fi
   else
      mkdir oldrsl_$it0".iter"$innerit0".stage1"
      mv rsl.* oldrsl_$it0".iter"$innerit0".stage1"
   fi
else
   mkdir oldrsl_$it0".iter"$innerit0_last".stage3"
   mv rsl.* oldrsl_$it0".iter"$innerit0_last".stage3"
fi

if [ $rand_stage -eq 3 ] && [ $rand_type -eq 3 ]; then
   #Check for presence of qhat vectors for this ensemble member (cv space)
   vectors=("qhat.e")
   vectors2=(".iter$innerit0")

   vcount=0
   for var in ${vectors[@]}
   do
      suffix=${vectors2[$vcount]}
      echo "Checking for $var$ii.*$suffix files"

      #Test for the presence of qhat for this ensemble member
      if [ $(ls ../vectors_$it0/$var$ii.*$suffix | wc -l) -eq 0 ]; then echo "ERROR: Missing $var$ii.*$suffix"; echo $((rand_stage*100+2)); exit $((rand_stage*100+2)); fi

      vcount=$((vcount+1))
   done
fi

#Handle external checkpoint and obs associated files that change each outer iteration
if [ $rand_stage -eq 1 ] && [ $rand_type -eq 3 ]; then
   if [ $CPDT -gt 0 ]; then
      rm wrf_checkpoint_d01_*
      if [ $WRF_MET -gt 0 ]; then
         rm xtraj_for_obs_d01_*
      fi
      if [ $WRF_CHEM -gt 0 ]; then
         rm SURFACE_Hx_y*
         rm AIRCRAFT_Hx_y*
      fi
   fi
else if ([ $rand_stage -eq 3 ] && [ $rand_type -eq 3 ]) || \
        ([ $rand_stage -eq 1 ] && [ $rand_type -eq 1 ]); then
   if [ $CPDT -gt 0 ]; then
      ln -sf $CWD_rel/wrf_checkpoint_d01_* ./
      if [ $WRF_MET -gt 0 ]; then
        ln -sf $CWD_rel/xtraj_for_obs_d01_* ./
      fi
   fi
   if [ $WRF_CHEM -gt 0 ]; then
      ln -sf $CWD_rel/SURFACE_Hx_y* ./
      ln -sf $CWD_rel/AIRCRAFT_Hx_y* ./
   fi
fi

#Link cvt file for current outer iteration
if [ $it -gt 1 ] && ([ $RIOT_RESTART -eq 0 ] || [ $RIOT_RESTART -eq 2 ]); then
   #Test for the presence of cvt
   if [ $(ls "$CWD"/cvt.it"$it0_last".* | wc -l) -eq 0 ]; then echo "ERROR: Missing cvt.*"; echo $((rand_stage*100+3)); exit $((rand_stage*100+3)); fi
   ln -sf $CWD_rel/cvt.* ./
fi


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

mpistring="$MPICALL $DEBUGSTR -np $npiens $RIOT_EXECUTABLE $BACKG_STRING"

# 2 - Potentially generic solution for Intel MPI implmentations
#      (may need earlier calls to mpdallexit and mpdboot)
#mpistring="$MPICALL -np $npiens -machinefile $(pwd)/hostlist $RIOT_EXECUTABLE $BACKG_STRING"

# 3 - Potentially generic solution for Open MPI implmentations
#      (may need earlier calls to mpdallexit and mpdboot)
#mpistring="$MPICALL -np $npiens --hostfile $(pwd)/hostlist $RIOT_EXECUTABLE $BACKG_STRING"

inithr0=$(($SECONDS / 3600))
if [ $inithr0 -lt 10 ]; then inithr0="0"$inithr0; fi
initmin0=$((($SECONDS / 60) % 60))
if [ $initmin0 -lt 10 ]; then initmin0="0"$initmin0; fi
initsec0=$(($SECONDS % 60))
if [ $initsec0 -lt 10 ]; then initsec0="0"$initsec0; fi

#Use these if $BACKGR_STRING does not place mpirun in background (no trailing "&"):
echo "JOB $iSAMP; $mpistring; $(date)"
eval "$mpistring"

##   #Use these if $BACKGR_STRING places mpirun in background (trailing "&"):
##   eval "$mpistring"; wait_pids+=($!)
##   echo "JOB $iSAMP; wait ${wait_pids[@]}"
##   wait "${wait_pids[@]}"

#NOTE: Redirected input ($BACKG_STRING) necessary to run in background and ensures clean log files

mpireturn=$?

hr0=$(($SECONDS / 3600))
if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
min0=$((($SECONDS / 60) % 60))
if [ $min0 -lt 10 ]; then min0="0"$min0; fi
sec0=$(($SECONDS % 60))
if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
echo "WRFVAR_MEMBER.$ii time: $hr0:$min0:$sec0;"$'\n'"INIT_TIME: $inithr0:$initmin0:$initsec0"$'\n'"w/ WRFDA return value: $mpireturn"

exit $mpireturn
