#!/bin/bash

#========================================================
#WRFVAR_ENSEMBLE.sh is used to run the ensemble for:
#stage 99 gradient preconditioning (rand_type==[1,6]) of WRFDA_RIOT.sh
#stage 1 (rand_type==[1,2,3,6]) of WRFDA_RIOT.sh
#stage 3 (rand_type==[1,3]) of WRFDA_RIOT.sh
#========================================================

#Probably should add a check for each of this I/O variables to make sure they are not empty
it="$1"; iSAMP0="$2"; NSAMP="$3"; NJOBS="$4"; rand_stage="$5"; innerit="$6"

it0=$it
if [ $it -lt 10 ]; then it0="0"$it0; fi

it0_last=$((it-1))
if [ $it0_last -lt 10 ]; then it0_last="0"$it0_last; fi

innerit0=$innerit
if [ $innerit -lt 10 ]; then innerit0=0$innerit0; fi
if [ $innerit -lt 100 ]; then innerit0=0$innerit0; fi
if [ $innerit -lt 1000 ]; then innerit0=0$innerit0; fi

proc_i=0
for (( iSAMP = $iSAMP0 ; iSAMP <= $NJOBS ; iSAMP++))
do
   echo "Starting job number $iSAMP for $NSAMP samples"
   ii=$iSAMP
   if [ $iSAMP -lt 10 ]; then ii=0$ii; fi
   if [ $iSAMP -lt 100 ]; then ii=0$ii; fi
   if [ $iSAMP -lt 1000 ]; then ii=0$ii; fi

   cd ../run.$ii

#   if [ $rand_stage -eq 1 ] && [ $iSAMP -le $NSAMP ] && \
#      [ "$GLOBAL_OMEGA" == "true" ] && \
#      ([ $it -eq 1 ] || [ $GRAD_PRECON -eq 0 ]); then
#      #Distribute omega vectors
#      #Test for the presence of each vector type
#      ls ../vectors_$it0/omega.e$ii.p0000
#      dummy=`ls ../vectors_$it0/omega.e$ii.p0000 | wc -l`
#      if [ $dummy -lt 1 ]; then
#         echo "ERROR: Missing global omega.e$ii.p0000"
#         echo $((rand_stage*100+1)); exit $((rand_stage*100+1))
#      fi
#      mv -v ../vectors_$it0/omega.e$ii* ./
#   fi

   if [ $innerit -eq 1 ]; then
      if [ $rand_stage -eq 1 ]; then
         if [ $it -gt 1 ] && \
            ([ $RIOT_RESTART -eq 0 ] || [ $RIOT_RESTART -eq 2 ]); then
            mkdir oldrsl_$it0_last
            mv -v rsl.* oldrsl_$it0_last
         fi
      else
         mkdir oldrsl_$it0".iter"$innerit0".stage1"
         mv -v rsl.* oldrsl_$it0".iter"$innerit0".stage1"
      fi
   else
      mkdir oldrsl_$it0".iter"$innerit0".stage3"
      mv -v rsl.* oldrsl_$it0".iter"$innerit0".stage3"
   fi


   if [ $rand_stage -eq 3 ]; then
      vectors2=("" "")
#      #Gather qhat_obs vectors for this ensemble (obs space)
#      if [ $rand_type -eq 1 ]; then vectors=("qhat_obs.e"); fi

      #Gather qhat vectors for this ensemble member (cv space)
      if [ $rand_type -eq 3 ]; then 
         vectors=("qhat.e")
         vectors2=(".iter$innerit0")
      fi

      vcount=0
      for var in ${vectors[@]}
      do
         echo "Working on $var files"
         suffix=${vectors2[$vcount]}

         #Test for the presence of qhat for this ensemble member
         if [ $(ls ../vectors_$it0/$var$ii.*$suffix | wc -l) -eq 0 ]; then echo "ERROR: Missing $var$ii.*$suffix"; echo $((rand_stage*100+2)); exit $((rand_stage*100+2)); fi

#         if [ $(ls $CWD_rel/$var$ii.*$suffix | wc -l) -eq 0 ]; then echo "ERROR: Missing $var$ii.*"; echo $((rand_stage*100+2)); exit $((rand_stage*100+2)); fi
#         mv -v $CWD_rel/$var$ii.*$suffix ../vectors_$it0
#
#         ln -sfv ../vectors_$it0/$var$ii.*$suffix ./

         vcount=$((vcount+1))
      done

#      if [ $innerit -gt 1 ]; then
         innerit_last=$((innerit-1))
         innerit_last0=$innerit_last
         if [ $innerit_last -lt 10 ]; then innerit0_last=0$innerit_last0; fi
         if [ $innerit_last -lt 100 ]; then innerit_last0=0$innerit_last0; fi
         if [ $innerit_last -lt 1000 ]; then innerit_last0=0$innerit_last0; fi

#      fi
   fi


   if [ $rand_stage -eq 1 ] || [ $rand_stage -eq 99 ]; then
      if [ $CPDT -gt 0 ]; then
         ln -sfv $CWD_rel/wrf_checkpoint_d01_* ./
         if [ $WRF_MET -gt 0 ]; then
           ln -sfv $CWD_rel/xtraj_for_obs_d01_* ./
         fi
      fi
      if [ $WRF_CHEM -gt 0 ]; then
         ln -sfv $CWD_rel/SURFACE_Hx_y* ./
         ln -sfv $CWD_rel/AIRCRAFT_Hx_y* ./
      fi

      if [ $it -gt 1 ] && ([ $RIOT_RESTART -eq 0 ] || [ $RIOT_RESTART -eq 2 ]); then
         #Test for the presence of cvt
         if [ $(ls "$CWD"/cvt.it"$it0_last".* | wc -l) -eq 0 ]; then echo "ERROR: Missing cvt.*"; echo $((rand_stage*100+3)); exit $((rand_stage*100+3)); fi
         ln -sfv $CWD_rel/cvt.* ./
      fi
   fi

   cp $CWD/namelist.input ./
   ex -c :"/ensmember" +:":s/=.*/=$iSAMP,/" +:wq namelist.input

#SHOULD SEPERATE INTO TWO LOOPS HERE TO CLOSE EFFICIENCY GAP (IDLE PROCESSES AFTER $RIOT_EXECUTABLE CALL)
   ## Assign the processes to the hostlist for the current ensemble member
   # Reverse consecutive process placement (tail)
   if [ $rand_stage -eq 99 ]; then
      if [ "$GLOBAL_OPT" == "true" ]; then
         npiens=$nproc_global
         #Could make this faster (more memory) by distributing across more nodes 
         # - currently chooses first $npiens processors in $PBS_NODEFILE
      else
         npiens=$nproc_local_grad
      fi
   else
      npiens=$nproc_local
      if [ $NSAMP -lt $NJOBS ] && [ $iSAMP -eq $NJOBS ]; then
         npiens=$nproc_local_grad
      fi
   fi

#DEBUG
#   echo "proc_f=$proc_f$"
#   echo "proc_f=$proc_i$"
#   echo "npiens=$npiens$"
#   echo "NUMPROC=$NUMPROC"
#   echo "PBSNODE0=$PBSNODE0"
#DEBUG
   proc_f=$((proc_i+npiens-1))
   tail -$((NUMPROC-proc_i)) $PBSNODE0 | head -$npiens > hostlist
   proc_i=$((proc_i+npiens))

   # Multiply A * w_i [iSAMP<=NSAMP] or A * (Hx - y) [iSAMP==NSAMP+1]
   # Note: The implementation through mpirun or mpiexec is 
   #       unique for your cluster and MPI implementation

   # DEFAULT - WORKS FOR:
   #  1a NASA Pleiades for SGI MPT
   #  1b NOAA Theia for Intel MPI
   #  *Please add more as verified or other solutions developed
   echo "export PBS_NODEFILE=$(pwd)/hostlist"
   export PBS_NODEFILE=$(pwd)/hostlist #.$ii
   #Redirected input necessary to run in background, output for clean log files
   mpistring="$MPICALL $DEBUGSTR -np $npiens $RIOT_EXECUTABLE $BACKG_STRING"

   # 2 - Generic solution for Intel MPI implmentations
   #      (may need earlier calls to mpdallexit and mpdboot)
   #mpistring="$MPICALL -np $npiens -machinefile $(pwd)/hostlist $RIOT_EXECUTABLE $BACKG_STRING"

   # 3 - Generic solution for Open MPI implmentations
   #      (may need earlier calls to mpdallexit and mpdboot)
   #mpistring="$MPICALL -np $npiens --hostfile $(pwd)/hostlist $RIOT_EXECUTABLE $BACKG_STRING"

#   if [ $rand_stage -eq 1 ] && [ $rand_type -eq 1 ] && \
#      [ $iSAMP -eq $NJOBS ] && [ $NSAMP -lt $NJOBS ]; then
#      eval "$mpistring"
#      echo "$mpistring"
#   else
      eval "$mpistring"; wait_pids+=($!)
      echo "$mpistring"
      echo "PID = "$!
#   fi

### Use these if eval above doesn't work
#         if [ $rand_type -eq 1 ] && [ $iSAMP -eq $NJOBS ]; then
#            $MPICALL $DEBUGSTR -np $npiens $RIOT_EXECUTABLE $BACKG_STRING
#
#            echo "PID = "$!
#            echo "$mpistring"
#         else
#            $MPICALL $DEBUGSTR -np $npiens $RIOT_EXECUTABLE $BACKG_STRING wait_pids+=($!)
#
#            echo "PID = "$!
#            echo "$mpistring"
#         fi

done

#WAIT for all ensembles to finish
echo "wait ${wait_pids[@]}"
wait "${wait_pids[@]}"

#Reset PBS_NODEFILE
export PBS_NODEFILE=$PBSNODE0

