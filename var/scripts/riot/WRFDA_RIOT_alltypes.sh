#!/bin/bash


SECONDS=0
iteration_sec=0
cleanup_sec=0
echo ""
echo "********************************************************"
echo ""
echo "         WRFDA-RIOT script with randomised SVD"
echo ""
echo "********************************************************"
echo ""

#=======================================================================================
# INSTRUCTIONS:
#=======================================================================================
# Place this script and RIOT_settings.sh in the same directory where you would normally 
# execute da_wrfvar.exe. Ensure that you are able to successfully run da_wrfvar.exe for 
# a 4D-Var inversion before using RIOT, including generation of cost_fn and 
# wrfvar_output files.
# The same namelist.input file used for standard 4D-Var will be used for RIOT 4D-Var,
#   with additional modifications automated by this script. You must have access to the
#   parent directory within which this working directory is located.  This script will
#   generate directories such as ../run.0001, ../run.0002, etc. that are essential.
# Fill in the environment variables in RIOT_settings.sh and the User Options section. 
#   - nout_RIOT takes the place of max_ext_its in namelist.input as the # outer iterations
#   - ntmax(1) in namelist.input is replaced by SVDN
#   - ntmax(2:nout_RIOT) should be set in namelist.input
#   - Follow the instructions for individual options below
# 
# Disclaimer: This script has been tested on a very small number of platforms. Each
#  platform will have a unique combination of MPI/MPT implementations and job queue 
#  management tools (e.g., PBS).  Specific solutions may be required for your 
#  environment. For questions/comments, please contact:
#
#  JJ Guerrette 
#  jonathan.guerrette@noaa.gov
#  jonathan.guerrette@colorado.edu
#
#=======================================================================================


if [ -z $MPICALL ]; then #Could be defined externally
   ##Check the manual for your compute system in order to choose one of the following:
   #MPICALL=mpiexec #(e.g., NASA Pleiades)
   MPICALL=mpirun #(e.g., NOAA Theia)

   # Note: the script performance has been verified on the platforms in parentheses;
   #       please add more as they are confirmed
fi
echo "MPI Calling Wrapper = "
echo "==> "$MPICALL

if [ -z "$WRF_CHEM" ]; then #Could be defined externally
   # Set to 0 (default) to not conduct CHEM DA
   # 1 to conduct chemical emission inversion (requires WRFDA-Chem)
   export WRF_CHEM=0
fi
if [ -z "$WRF_MET" ]; then #Could be defined externally
   # Set to 0 (default) to not conduct MET DA
   # 1 to conduct MET DA
   export WRF_MET=1
fi
if [ $WRF_CHEM -eq 0 ] && [ $WRF_MET -eq 0 ]; then
   echo "ERROR: Either WRF_CHEM or WRF_MET must be set to 1"
   echo 1; exit 1
fi
#=======================================================================================
# Begin User Options
#=======================================================================================
#####################################################################################
#Modify RIOT_settings.sh to control RIOT behavior

#Read in baseline RIOT settings
if [ -z $RIOT_SETTINGS_CALLED ]; then
   . ./RIOT_settings.sh
fi

#Manually set the maximum number of processes per job 
# - limited by WRF patch overlap
# - depends on domain (nx x ny) and PPN
# - critical for speeding up single-job (serial) portions of RIOT 
#    (e.g., GRAD_PRECON=1, STAGE 2 of RSVD5.6 and RSVD5.1, and STAGE 4 of RSVD5.1)
# - A good rule of thumb is (nx/10 * ny/10) <= NPpJMAX << (nx/5 * ny/5)
# - Requirement: NPpJMAX <= NUMPROC (see below)
if [ -z $NPpJMAX ]; then
   export NPpJMAX=64
fi

#Manually turn on/off gathering of subprocedure times
export SUBTIMING=1 #0 or 1
#####################################################################################

#Additional local settings
OUTDIR="ALL_RIOT_OUTPUT"

#MPI options
DEBUGSTR=
#DEBUGSTR="-show"
#DEBUGSTR="-verbose"

export MPI_VERBOSE=1
export BACKG_STRING="< /dev/null &> out.run &"
export EXECUTABLE="./da_wrfvar.exe"

#------------------------------------
CLEANUP=1 #Takes extra time or space:
# 0 - leave temp files as-is
# 1 - store temp files in tars
# >1 - remove temp files
#------------------------------------

#All of these must be set to 1 to perform RIOT - toggle for debugging
STAGE0=1
STAGE1=1 #Set to 1 to perform parallel A * OMEGA calculations
STAGE2=1 #Set to 1 to perform final SVD and increment CVs

#=======================================================================================
# End User Options
#=======================================================================================

echo "=================================================="
echo " Critical RIOT Options"
echo "=================================================="
#Retrieve ntmax array from namelist - set NENS for it>1
ntmax_all=`grep ntmax namelist.input`
ntmax_all=${ntmax_all#*=}
IFS=',' read -ra ntmax_array <<< "$ntmax_all"

echo "(1) rand_type=$rand_type"
echo " * valid rand_type options: "
echo "   + 6-RSVD5.6"
echo "   + 1-RSVD5.1"
echo "   + 2-Full Hessian(Requires SVDN=Nobs)"
echo "   + 10-B cov. only, parallel (get leading eigenmodes)"
echo "   + 11-B cov. only, serial (get leading eigenmodes)"
echo " * only rand_type=6, 10, 11 are functional in WRFDA for non-CHEM variables"
if [ $rand_type -ne 1 ] && [ $rand_type -ne 6 ] && [ $rand_type -ne 2 ] && [ $rand_type -ne 10 ] && [ $rand_type -ne 11 ]; then echo "ERROR: unknown rand_type=$rand_type"; echo 1; exit 1; fi
echo ""
echo "(2) RIOT_RESTART=$RIOT_RESTART"
echo ""
echo "(3) ALT_START=$ALT_START"
echo ""
echo "(4) MAX processors per job: NPpJMAX=$NPpJMAX"
echo ""
#Could use these two lines instead:
#NOUTER=`grep max_ext_its namelist.input`
#NOUTER=${NOUTER//[!0-9]/}
NOUTER=$nout_RIOT # overwrites NOUTER with environment variable
echo "(5) # of outer iterations: nout_RIOT=$NOUTER"
ex -c :"/max_ext_its" +:":s/=.*/=1,/" +:wq namelist.input

NENS=$SVDN
ntmax_array[0]=$SVDN
for it in $(seq 1 $nout_RIOT)
do
   echo "# ensembles in outer iteration $it: ntmax($it)=${ntmax_array[$((it-1))]}"
   ntmax=$ntmax${ntmax_array[$((it-1))]}","
done
ntmax=$ntmax${ntmax_array[$((nout_RIOT-1))]}","
ex -c :"/ntmax" +:":s/=.*/=$ntmax/" +:wq namelist.input

echo ""
echo "(7) LRA-LRU Adaptation: ADAPT_SVD=$ADAPT_SVD"
#if [ "$ADAPT_SVD" != "false" ] && [ "$ADAPT_SVD" != "true" ]; then
#   echo "ERROR: ADAPT_SVD must either be true or false"; echo 2; exit 2
#fi
echo ""
echo "(8) Preconditioning Option: GRAD_PRECON=$GRAD_PRECON"
valid_precon=0
for i in 1 2 3 4 12 13 14 15
do
   if [ "$GRAD_PRECON" -eq $i ]; then valid_precon=1; fi
done
if [ "$GRAD_PRECON" -lt 0 ] || [ $valid_precon -eq 0 ]; then
#if [ "$GRAD_PRECON" -lt 0 ] || [ "$GRAD_PRECON" -gt 4 ]; then
   echo "ERROR: GRAD_PRECON must be >= 0 and valid values are"
   echo " 0, 1, 2, 3, 4, 12, 13, 14, or 15"; echo 3; exit 3
fi
echo ""
if [ -z $GLOBAL_OPT ]; then
   GLOBAL_OPT="true" #default choice
fi
echo "(9) GLOBAL_OPT=$GLOBAL_OPT"
if [ "$GLOBAL_OPT" == "true" ]; then
   echo "For very large size(cv), global cv I/O requires"
   echo "1 - additional wall-time during WRFDA gather step"
   echo "2 - larger minimum memory on head node, for same step"
   echo "Consider these factors if da_wrfvar.exe hangs or quits unexpectedly during cv I/O"
   echo "Turn off global I/O with, e.g., export GLOBAL_OPT=\"false\""
   echo ""
else if [ "$GLOBAL_OPT" == "false" ]; then
   echo "Note: Global I/O is necessary to change number of processors per ensemble member between outer iterations"
else
   echo "ERROR: GLOBAL_OPT must either be true or false"; echo 4; exit 4
fi
fi

echo "================================================="
echo "Characerize the nodes and processors being used"
echo "================================================="
PPN=  #Use this to override the number of processors per node
if [ -z $PBS_NODEFILE ]; then
   echo "ERROR: Using SVD requires PBS_NODEFILE to be set to your host file."
   echo "ERROR: This tool requires a large number of processors (> 100)"
   echo 5; exit 5
fi
PBSNODE0=$PBS_NODEFILE
NUMNODES=`cat $PBSNODE0 | uniq | wc -l`
if [ $NUMNODES -lt $((SVDN+1)) ]; then
   echo "ERROR: NUMNODES must be set >= SVDN+1 (i.e., $((SVDN+1)))"
   echo 6; exit 6
fi

if [ -n "$PPN" ]; then
   echo "Using user defined cores per node: PPN=$PPN"
   NUMPROC=$((NUMNODES*PPN))
else
   NUMPROC=`cat $PBSNODE0 | wc -l`
   PPN=$((NUMPROC/NUMNODES))
   echo "Automatically determining cores per node: PPN=$PPN"
fi

if [ $NPpJMAX -gt $NUMPROC ]; then
   NPpJMAX=$NUMPROC
   export NPpJMAX=$NPpJMAX
#   echo "ERROR: NPpJMAX must be set <= total number of processors"
#   echo 7; exit 7
fi
#Use as many of the cores as possible for ensemble members
#NODES_per_ens=$((NUMNODES/NENS))
#NPpJ=$((NODES_per_ens*PPN))

NPpJ=$(($((NUMNODES/NENS))*PPN)) 
if [ $NPpJ -gt $NPpJMAX ]; then NPpJ=$(($((NPpJMAX/PPN))*PPN)); fi
NODES_per_ens=$((NPpJ/PPN*NENS))

#Use remaining cores for gradient member
NODES_grad=$((NUMNODES-NODES_per_ens))
#Fail-safe for when mod(NUMNODES,NENS)==0
if [ $NODES_grad -lt 1 ]; then
   #NODES_per_ens=$(($((NUMNODES-1))/NENS))
   #NPpJ=$((NODES_per_ens*PPN))
   NPpJ=$(($(($((NUMNODES-1))/NENS))*PPN)) 
   NODES_per_ens=$((NPpJ/PPN*NENS))
   NODES_grad=$((NUMNODES-NODES_per_ens))
fi
NPpJ_grad=$((NODES_grad*PPN))
if [ $NPpJ_grad -gt $NPpJMAX ]; then NPpJ_grad=$(($((NPpJMAX/PPN))*PPN)); fi
NODES_grad=$((NPpJ_grad/PPN))

if [ "$GLOBAL_OPT" == "false" ] && [ $NPpJ_grad -ne $NPpJ ]; then
    echo "ERROR: nproc_local and nproc_local_grad must be equal when GLOBAL_OPT==false"
    echo "NPpJ=$NPpJ"
    echo "NPpJ_grad=$NPpJ_grad"
    echo 8; exit 8
fi

echo "Cores per ensemble job: "$NPpJ
echo "Cores per non-global gradient job: "$NPpJ_grad
if [ $((NPpJ % PPN)) -ne 0 ] || [ $((NPpJ_grad % PPN)) -ne 0 ] || [ $NPpJ -eq 0 ] || [ $NPpJ_grad -eq 0 ]; then
    echo "ERROR: cores per ensemble must be equal to PPN times a positive integer multiplier"
    echo "NPpJ=$NPpJ"
    echo "NPpJ_grad=$NPpJ_grad"
    echo "PPN=$PPN"
    echo 9; exit 9
fi

# Set processor counts for: 
#(1) stage 0
nproc_max=$NUMPROC
if [ $nproc_max -gt $NPpJMAX ]; then nproc_max=$NPpJMAX; fi

#(2) ensemble members
nproc_local=$NPpJ

#(3) gradient member
nproc_local_grad=$NPpJ_grad

#(4) global gathering stages
if [ "$GLOBAL_OPT" == "false" ]; then
   #FOR NON-GLOBAL SIM
   nproc_global=$nproc_local
else if [ "$GLOBAL_OPT" == "true" ]; then
   #FOR GLOBAL ENABLED SIM
   nproc_global=$nproc_max
fi
fi

echo "LOCAL PROCS = $nproc_local"
echo "GRAD PROCS = $nproc_local_grad"
echo "GLOBAL PROCS = $nproc_global"

# Set number of parallel jobs
NJOBS=$NENS
if [ $rand_type -ne 10 ] && [ $rand_type -ne 11 ]; then
   NJOBS=$((NENS+1)) #Extra parallel ensemble for gradient member
fi


echo "=================================================================="
echo " Automatically initialize 4D-Var namelist.input settings for RIOT"
echo "=================================================================="

#------------------------------------------------------------------
# Set checkpoint interval and OSSE from environment variables
# - only meaningful for CHEM so far
# - set to 0 otherwise
# - if non-zero, namelist values should already be set
#------------------------------------------------------------------
CPDT=$checkpoint_interval 
quick_rand="true"
if [ -z $CPDT ]; then
   if [ $WRF_CHEM -eq 0 ]; then
      CPDT=0 #Or set to some other default for WRF_CHEM==0
   else
      #Thus far WRF_CHEM==1 requires CPDT>0
      CPDT=180 # (or manually, e.g., 5, 6, 10, 12, 15, 20, 30, 60, 180)
   fi
fi
if [ $CPDT -le 0 ]; then quick_rand="false"; fi
echo "checkpoint interval = "
echo "==> "$CPDT

#CHEM OSSE

OSSE=$OSSE_CHEM # Set OSSE (from env var)
if [ -z $OSSE ] || [ $WRF_CHEM -eq 0 ]; then
   OSSE=0 #Default if env var missing
fi
echo "OSSE = "
echo "==> "$OSSE


# Turn off lanczos
if (grep -q use_lanczos namelist.input); then
   ex -c :"/use_lanczos" +:":s/use_lanczos.*/use_lanczos=false,/" +:wq namelist.input
fi

# Add or modify RIOT namelist options
SVD_VARS=("rand_outer" "rand_minimise" "ensmember" "rand_stage" "use_randomblock" "rand_type" "quick_rand" "adapt_svd" "svd_p" "riot_precon" "read_omega" "use_global_cv_io" "prepend_rsvd_basis")
SVD_VALS=("1" "true" "0" "0" "true" "$rand_type" "$quick_rand" "$ADAPT_SVD" "$svd_p" "$GRAD_PRECON" "$GLOBAL_OMEGA" "$GLOBAL_OPT" "0")
ivar=0
for var in ${SVD_VARS[@]}
do
   if (grep -q $var namelist.input); then # check if NL option is present already
      ex -c :"/$var" +:":s/=.*/=${SVD_VALS[$ivar]},/" +:wq namelist.input
   else
      ex -c :"/wrfvar6" +:":s/\(wrfvar6\)\n/\1\r$var=${SVD_VALS[$ivar]},\r/" +:wq namelist.input
   fi
   ivar=$((ivar+1))
done
prepend_basis=" 0," #For now, do not prepend in first iteration
for it in $(seq 2 $nout_RIOT)
do
   prepend_basis=$prepend_basis" $prepend_rsvd_basis,"
done
ex -c :"/prepend_basis" +:":s/=.*/=$prepend_basis/" +:wq namelist.input


echo "=================================================="
echo " Setup ensemble directories"
echo "=================================================="
CWD=$(pwd)
echo "WORKING DIRECTORY: "$CWD
cd ../
PARENTDIR=$(pwd)
CWD_rel=${CWD#$PARENTDIR"/"} #Use for linking, more stable
CWD_rel="../"$CWD_rel
cd $CWD

dummy=$(ls ../run.*)
if [ $? -eq 0 ]; then  rm -r ../run.[0-9]*; fi

if [ $OSSE -eq 2 ]; then
#   DIRpert="LZ_N=100_no=1_LARGE"
   DIRpert="LZ_N=10_no=10_LARGE"
   fg_pert="/nobackupp8/jjguerr1/wrf/DA/"$DIRpert"/run/fg_pert"
   osse_obs="/nobackupp8/jjguerr1/wrf/DA/"$DIRpert"/run/AIRCRAFT_DATA_*.txt"
   osse_var="/nobackupp8/jjguerr1/wrf/DA/"$DIRpert"/run/AIRCRAFT_MODEL_VARIANCE_*.txt"

   ln -sf $fg_pert ./fg_pert
   ln -sf $osse_obs ./
   ln -sf $osse_var ./
fi

if [ $RIOT_RESTART ] && [ $RIOT_RESTART -eq 1 ]; then
#----------------------------------------------------------------
# Option to perform posterior error calculation w/o minimization
#----------------------------------------------------------------
# Either setup links here or externally
#   ln -sf $ALT_START ./fg
#   ln -sf $ALT_START ./wrfinput_d01

   itstart=2
   NOUTER=$((NOUTER+1))
   ex -c :"/rand_minimise" +:":s/=.*/=false,/" +:wq namelist.input
else if [ $RIOT_RESTART ] && [ $RIOT_RESTART -eq 2 ]; then
#------------------------------------------------------------------
# Option to restart RIOT minimization from previous outer iteration
#------------------------------------------------------------------
# Either setup links here or externally
#   ln -sf $ALT_START ./fg
#   ln -sf $ALT_START ./wrfinput_d01

   itstart=$ALT_it1
   if [ -z $itstart ]; then
      itstart=2
   fi
   ii=$((itstart-1))
   if [ $ii -lt 10 ]; then ii=0$ii; fi
   ln -sf $ALT_CVT ./cvt.it"$ii".p0000
   ln -sf $ALT_XHAT ./xhat.it"$ii".p0000

   if [ $GRAD_PRECON -gt 0 ]; then
      cd ../
      ln -sf $ALT_hess_dir/hessian_eigenpairs.it* ./
#      ln -sf $ALT_hess ./hessian_eigenpairs.it"$ii".0000
      cd $CWD
   fi
else
   itstart=1 #default
fi
fi

#Populate ensemble directories
if [ $rand_type -ne 11 ]; then
   if [ $(ls "../run.*" | wc -l) -gt 0 ]; then  rm -r ../run.[0-9]*; fi
   if [ $WRF_CHEM -gt 0 ]; then rm *Hx*; fi

   for (( iENS = 1 ; iENS <= $NJOBS ; iENS++))
   do
      ii=$iENS
      if [ $iENS -lt 10 ]; then ii=0$ii; fi
      if [ $iENS -lt 100 ]; then ii=0$ii; fi
      if [ $iENS -lt 1000 ]; then ii=0$ii; fi

      mkdir ../run.$ii
      cd ../run.$ii

      # link all files necessary to run da_wrfvar.exe
      ln -s $CWD_rel/* ./
      rm namelist.input
      rm rsl.out.*
      rm rsl.error.*
   done
fi
cd $CWD

hr0=$(($SECONDS / 3600))
if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
min0=$((($SECONDS / 60) % 60))
if [ $min0 -lt 10 ]; then min0="0"$min0; fi
sec0=$(($SECONDS % 60))
if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
echo "RIOT initialization time: $hr0:$min0:$sec0"

#-----------------------------------------
#Export variables used in WRFVAR_ENSEMBLE
#-----------------------------------------
export CWD
export CWD_rel
export rand_type
export GLOBAL_OMEGA
export GLOBAL_OPT
export GRAD_PRECON
export CPDT
export WRF_CHEM
export RIOT_RESTART
export nproc_local
export nproc_local_grad
export nproc_global
export NUMPROC
export PBSNODE0
export MPICALL
export DEBUGSTR
export EXECUTABLE
export BACKG_STRING

if [ $rand_type -eq 10 ]; then
   itstart=1
   NOUTER=1
fi

## BEGIN Outer Loop
for (( it = $itstart ; it <= $NOUTER ; it++))
do
   SECONDS=0
   echo ""
   echo ""
   echo "===================================================="
   echo "===================================================="
   echo "==                                                =="
   echo "==  Starting outer loop iteration $it of $NOUTER  =="
   echo "==                                                =="
   echo "===================================================="
   echo "===================================================="
   echo ""
   echo ""
   ex -c :"/rand_outer" +:":s/=.*/=$it,/" +:wq namelist.input
   it0=$it
   if [ $it -lt 10 ]; then it0="0"$it0; fi

   cd $CWD

   if [ $(ls -d "../vectors_$it0" | wc -l) -gt 0 ]; then rm -rv ../vectors_$it0; fi
   mkdir -v ../vectors_$it0
   rm omega.e*.p*
   rm yhat.e*.p*

   it0_last=$((it-1))

   NENS=${ntmax_array[$((it-1))]}

   #Hard constraint that NENS must be <= SVDN (necessary?)
   if [ $NENS -gt $SVDN ]; then
      echo "WARNING SVDN=$SVDN, ntmax(it)=$NENS, setting NENS=$SVDN"
      NENS=$SVDN
      ntmax_array[$((it-1))]=$NENS
      ntmax=
      for it in $(seq 1 $nout_RIOT)
      do
         echo "# ensembles in outer iteration $it: ntmax($it)=${ntmax_array[$((it-1))]}"
         ntmax=$ntmax${ntmax_array[$((it-1))]}","
      done
      ntmax=$ntmax${ntmax_array[$((nout_RIOT-1))]}","
      ex -c :"/ntmax" +:":s/=.*/=$ntmax/" +:wq namelist.input
   fi
   echo ""
   echo "Using ensemble size $NENS in iteration $it, according to ntmax in namelist.input."
   echo ""

   if [ $it0_last -lt 10 ]; then it0_last="0"$it0_last; fi
   if [ $STAGE0 -gt 0 ]; then
#====================================================================================
#====================================================================================
      echo "============="
      echo "SVD STAGE 0"
      echo "============="
      # This stage only necessary to generate checkpoint or OSSE files (CHEM only)

      ex -c :"/rand_stage" +:":s/=.*/=0,/" +:wq namelist.input
      if [ $OSSE -gt 0 ] && [ $it -eq 1 ]; then
         echo "Generating/linking OSSE files..."
         if [ $OSSE -eq 1 ]; then
            ex -c :"/osse_chem" +:":s/=.*,/=false,/g" +:wq namelist.input
            ex -c :"/init_osse_chem" +:":s/=.*,/=true,/g" +:wq namelist.input
            mpistring="$MPICALL $DEBUGSTR -np $nproc_max $EXECUTABLE"
            #COULD MAKE THIS FASTER (MORE MEMORY) by distributing across more nodes (currently chooses first $npiens processors in $PBS_NODEFILE)

            echo "$mpistring"
            eval "$mpistring"

            mkdir rsl_init_osse
            cp rsl.* rsl_init_osse

            mv fgc fgc_0
         fi
         ln -sf ./fg_pert ./fg
         ln -sf ./fg_pert ./wrfinput_d01

         ex -c :"/osse_chem" +:":s/=.*,/=true,/g" +:wq namelist.input
         ex -c :"/init_osse_chem" +:":s/=.*,/=false,/g" +:wq namelist.input
      fi

      if [ $CPDT -gt 0 ] || [ "$GLOBAL_OMEGA" == "true"  ] ; then
         echo "Generating checkpoint files and/or global omega..."

         if [ "$GLOBAL_OMEGA" == "true" ] && [ "$GLOBAL_OPT" == "false" ]; then
            npiens=$nproc_local
         else
            npiens=$nproc_global
               #COULD MAKE THIS FASTER (MORE MEMORY) by distributing across more nodes (currently chooses first $npiens processors in $PBS_NODEFILE)
         fi

         mpistring="$MPICALL $DEBUGSTR -np $npiens $EXECUTABLE"
#         mpistring="$MPICALL $DEBUGSTR -np $nproc_max $EXECUTABLE"
         echo "$mpistring"
         eval "$mpistring"

         if [ $SUBTIMING -eq 1 ]; then grep da_end_timing rsl.out.0000 > ../bench_time_stage0.it$it0; fi
      fi
   fi

   if [ $STAGE1 -le 0 ] || [ $rand_type -eq 11 ]; then
      echo "EXIT: STAGE1 > 0 required for multiple outer iterations"; echo 10; exit 10;
   fi

   if [ "$GLOBAL_OMEGA" == "true" ]; then
      if [ $it -eq 1 ] || [ $GRAD_PRECON -eq 0 ]; then
         echo "Distributing global omega files"

         #Test for the presence of each vector type
         cd $CWD
         ls omega.e*.p0000
         dummy=`ls omega.e*.p0000 | wc -l`
         echo "$dummy present of $NENS omega files"
         if [ $dummy -ne $NENS ]; then 
            echo "ERROR: Missing or extra omega.e*.p0000"
            echo 18; exit 18
         fi
         mv -v omega.e*.p* ../vectors_$it0
      else
         rm omega.e*.p*
      fi
   fi


#====================================================================================
#====================================================================================
   echo "======================================================="
   echo "SVD STAGE 1: Multiply A * OMEGA, and calculate GRAD(J)"
   echo "======================================================="
   ex -c :"/rand_stage" +:":s/=.*/=1,/" +:wq namelist.input

   # Check for presence of checkpoint and obs output files (CHEM only)
   if [ $CPDT -gt 0 ]; then
      if [ $(ls "$CWD"/wrf_checkpoint_d01* | wc -l) -eq 0 ]; then echo "ERROR: Missing checkpoint files"; echo 11; exit 11; fi
   fi
   if [ $WRF_CHEM -gt 0 ]; then
      if [ $(ls "$CWD"/SURFACE_Hx_y* | wc -l) -eq 0 ]; then echo "ERROR: Missing SURFACE_Hx_y files"; echo 12; exit 12; fi

      if [ $(ls "$CWD"/AIRCRAFT_Hx_y* | wc -l) -eq 0 ]; then echo "ERROR: Missing AIRCRAFT_Hx_y files"; echo 13; exit 13; fi
# Something like this would allow storing checkpoint files in a temporary directory (hard disk, not memory due to size)
#         if [ $it -gt $itstart ]; then
#            pbsdsh -- rm $TMPDIR/wrf_checkpoint_d01*
#         fi
#
#         #These cp's take extra time, which could be avoided if simultaneous access is allowed.
#         pbsdsh -- cp $CWD/wrf_checkpoint_d01* $TMPDIR
#         pbsdsh -- cp $CWD/SURFACE_Hx_y* $TMPDIR
#         pbsdsh -- cp $CWD/AIRCRAFT_Hx_y* $TMPDIR
   fi

   if [ $it -gt $itstart ] || ([ $it -gt 1 ] && [ $RIOT_RESTART -eq 2 ]); then
      NENS_PREV=0
      for itprev in $(seq 1 $((it-1)))
      do
#         NENS_last=${ntmax_array[$((itprev-1))]}
         NENS_PREV=$((NENS_PREV+${ntmax_array[$((itprev-1))]}))
         echo "NENS_PREV=$NENS_PREV"
      done

      # This test avoids extra wall time when GRAD_PRECON=[1,3] and NENS in 
      # current iteration is larger than number of potential preconditioning vectors
      dependent_grad=0
      if [ $NENS_PREV -gt $NENS ] && \
        ([ $GRAD_PRECON -eq 1 ] || [ $GRAD_PRECON -eq 3 ]); then dependent_grad=1; fi
      if [ $GRAD_PRECON -eq 2 ] || [ $GRAD_PRECON -eq 4 ]; then dependent_grad=1; fi

      NJOBS=$((NENS+1))
      ii_grad=$NJOBS
      if [ $NJOBS -lt 10 ]; then ii_grad=0$ii_grad; fi
      if [ $NJOBS -lt 100 ]; then ii_grad=0$ii_grad; fi
      if [ $NJOBS -lt 1000 ]; then ii_grad=0$ii_grad; fi

      if [ $rand_type -eq 10 ]; then
         NJOBS=$NENS
      fi

      if [ $GRAD_PRECON -gt 0 ] && ([ $RIOT_RESTART -eq 0 ] || [ $RIOT_RESTART -eq 2 ]); then
#====================================================================================
         echo ""
         echo "------------------------------------------------------"
         echo "Setting up RIOT preconditioning for outer loop = $it"
         echo "------------------------------------------------------"
         echo ""
         ex -c :"/read_omega" +:":s/=.*/=true,/" +:wq namelist.input

         #Calculate gradient and preconditioned OMEGA vectors first
         ii=$ii_grad

#         if ([ $NENS_PREV -gt $NENS ] \
#            && ([ $GRAD_PRECON -eq 1 ] \
#             || [ $GRAD_PRECON -eq 3 ])) \
#             || [ $GRAD_PRECON -eq 2 ] \
#             || [ $GRAD_PRECON -eq 4 ] ; then
#            dependent_grad=1; else dependent_grad=0; fi

         if [ $dependent_grad -eq 1 ]; then
            #Calculate gradient and preconditioned OMEGA vectors before ensemble
            echo "Calulating gradient-eignevector dot products, "
            echo "and selecting $NENS preconditioned perturbation "
            echo "members (omega_i=eignevector_i) with largest values."

            ./WRFVAR_ENSEMBLE.sh $it $NJOBS $NENS $NJOBS 3

            diri=../run.$ii

            NJOBS=$NENS
         else
            #GRAD_PRECON=12,13,14,15 --> EXTRACT omega_precon.e*.p* only
            diri=$CWD
         fi

         echo ""
         echo "Transfer preconditioned OMEGA vectors before STAGE 1..."
         for (( iENS = 1 ; iENS <= $NENS ; iENS++))
         do
            jj=$iENS
            if [ $iENS -lt 10 ]; then jj=0$jj; fi
            if [ $iENS -lt 100 ]; then jj=0$jj; fi
            if [ $iENS -lt 1000 ]; then jj=0$jj; fi

            cd ../run.$jj
            ls $diri/omega_precon.e$jj.p*
            if [ $? -ne 0 ]; then echo "ERROR: Missing run.$ii/omega.e$jj.p*"; echo 15; exit 15; fi

            for file in $diri/omega_precon.e$jj.p*
            do
               myfile=${file/_precon}
               myfile=${myfile#$diri/}
               mv -v $file ./$myfile
            done
         done
      fi
#====================================================================================

      #----------------------------------------------------
      # Redistribute cores among ensembles as NENS changes
      #----------------------------------------------------
      # When max(ntmax_array(iteration>=it)) ~= max(ntmax_array(iteration<it)), 
      # or ntmax_array(it) ~= ntmax_array(it-1), there are several
      # options for improved efficiency or performance:

      ##--------------------------------------------------------------------------
      ## 1 - Release nodes that would go unused (narrow the job width on the fly)
      ##  --> Possible on NAS Pleiades with pbs_release_nodes...still need to script
      #
      ##--------------------------------------------------------------------------

      ##--------------------------------------------------------------------------
      ##2 - Adjust cores per job to fit NUMPROC as NENS changes
      if [ "$GLOBAL_OPT" == "true" ]; then # && [ $NENS -ne $NENS_last ]; then
         if [ $dependent_grad -eq 1 ]; then
            NODES_per_ens=$NUMNODES
         else
            NODES_per_ens=$((NUMNODES-NODES_grad))
         fi
         NPpJ=$(($((NODES_per_ens/NENS))*PPN)) #Use remaining cores for gradient member
         if [ $((NPpJ/PPN)) -eq 0 ]; then
            echo "ERROR: cores per ensemble less than cores per node"
            echo "$NUMPROC, $NODES_per_ens, $PPN, $NPpJ"
            echo 16; exit 16
         fi
         if [ $NPpJ -gt $NPpJMAX ]; then NPpJ=$(($((NPpJMAX/PPN))*PPN)); fi

         echo "Cores per ensemble job: "$NPpJ
         nproc_local=$NPpJ

         dummy=$(($((nproc_local*NENS))+$(($((NJOBS-NENS))*nproc_local_grad))))
         if [ $dummy -gt $NUMPROC ]; then
            echo "ERROR: Too many processors requested in it=$it, NUMPROC_request=$dummy, NUMPROC_avail=$NUMPROC"; echo 17; exit 17
         fi
      fi
      ##--------------------------------------------------------------------------

      ##--------------------------------------------------------------------------
      ##3 - Some combination of these two --> custom option for particular applications
      #
      ##--------------------------------------------------------------------------

      #2 and #3 are possible only with global CV I/O in WRFDA (GLOBAL_OPT="true")
   fi

   #==============================================
   # INITIATE STAGE1 ENSEMBLE FOR ALL SVD TYPES
   #==============================================
   ./WRFVAR_ENSEMBLE.sh $it 1 $NENS $NJOBS 3

   #WAIT for all ensembles to finish
   echo "wait ${wait_pids[@]}"

   wait "${wait_pids[@]}"
#      wait

   #Reset PBS_NODEFILE
   export PBS_NODEFILE=$PBSNODE0

   if [ $STAGE2 -le 0 ]; then
      echo "EXIT: STAGE2 > 0 required for multiple outer iterations"; echo 21; exit 21;
   fi
#===================================================================================
   if [ $rand_type -eq 1 ]; then
      echo "=================================================="
      echo "SVD STAGE 2: Generating observation space basis, Q"
      echo "=================================================="
   fi
   if [ $rand_type -eq 2 ]; then
      echo "======================================================================"
      echo "SVD STAGE 2: Perform Eigen Decomp + Increment CVs using direct Hessian"
      echo "======================================================================"
   fi
   if [ $rand_type -eq 6 ]; then
      echo "================================================="
      echo "SVD STAGE 2: Perform Eigen Decomp + Increment CVs"
      echo "================================================="
   fi
   if [ $rand_type -eq 10 ]; then
      echo "================================================="
      echo "SVD STAGE 2: Perform Eigen Decomp of B Matrix"
      echo "================================================="
   fi
#===================================================================================

   cd $CWD

   #Collect relevant vectors into single directory
   if [ $rand_type -eq 1 ]; then
      #Gather yhat_obs vectors (obs space)
      vectors=("yhat_obs.e*.p")
      nvec=("$NENS")
   fi
   if [ $rand_type -eq 2 ]; then
      #Gather ahat and ghat vectors (cv space)
      vectors=("ahat.e*.p" "ghat.p")
      nvec=("$NENS" "1")
   fi
   if [ $rand_type -eq 6 ]; then
      #Gather omega, yhat, and ghat vectors (cv space)
      vectors=("omega.e*.p" "yhat.e*.p" "ghat.p")
      nvec=("$NENS" "$NENS" "1")
   fi
   if [ $rand_type -eq 10 ]; then
      #Gather omega, yhat, and ghat vectors (cv space)
      vectors=("omega.e*.p" "yhat.e*.p")
      nvec=("$NENS" "$NENS" "1")
   fi

   vcount=0
   for var in ${vectors[@]}
   do
      echo "Working on $var files"
      #Test for the presence of each vector type
      ls ../run.*/$var"0000"
      dummy=`ls ../run.*/$var"0000" | wc -l`
      echo "$dummy present of ${nvec[$vcount]} $var files"
      if [ $dummy -ne ${nvec[$vcount]} ]; then 
         echo "ERROR: Missing or extra $var""0000"
         echo 22; exit 22
      fi
      mv -v ../run.*/$var* ../vectors_$it0
      ln -sf ../vectors_$it0/$var* ./
      vcount=$((vcount+1))
   done

   ex -c :"/rand_stage" +:":s/=.*/=2,/" +:wq namelist.input

   if [ $rand_type -eq 1 ]; then
      npiens=$nproc_local  #can remove this when yhat_obs (and qhat_obs) has global I/O
   else
      npiens=$nproc_global
   fi

   #Perform Eigen Decomp + Calculate Increment and Analysis
   mpistring="$MPICALL $DEBUGSTR -np $npiens $EXECUTABLE"
   #COULD MAKE THIS FASTER (MORE MEMORY) by distributing across more nodes (currently chooses first $npiens processors in $PBS_NODEFILE)

   echo "$mpistring"
   eval "$mpistring"

   if [ $rand_type -eq 6 ] && [ $SUBTIMING -eq 1 ]; then
      grep da_end_timing ../run.*/rsl.out.0000 > ../bench_time_hess-vec.it$it0
      grep da_end_timing rsl.out.0000 > ../bench_time_finalize-riot.it$it0
   fi

   if [ $rand_type -eq 1 ]; then
#===================================================================
      echo "-------------------------------------------------"
      echo "SVD STAGE 3: Evaluating KMAT^T = A^T * Q"
      echo "-------------------------------------------------"
#===================================================================
      ex -c :"/rand_stage" +:":s/=.*/=3,/" +:wq namelist.input

      ./WRFVAR_ENSEMBLE.sh $it 1 $NENS $NENS 3

      #Reset PBS_NODEFILE
      export PBS_NODEFILE=$PBSNODE0


#===================================================================
      echo "-------------------------------------------------"
      echo "SVD STAGE 4: Performing SVD of KMAT = Q^T * A"
      echo "             + Increment CVs"
      echo "-------------------------------------------------"
#===================================================================

      cd $CWD
      if [ $SUBTIMING -eq 1 ]; then grep da_end_timing ../run.*/rsl.out.0000 > ../bench_time_hess-vecA.it$it0; fi

      #Gather bhat and ghat vectors (cv space)
      vectors=("bhat.e*.p" "ghat.p")
      nvec=("$NENS" "1")
      vcount=0
      for var in ${vectors[@]}
      do
         echo "Working on $var files"
         #Test for the presence of each vector type
         ls ../run.*/$var"0000"
         dummy=`ls ../run.*/$var"0000" | wc -l`
         echo "$dummy present of ${nvec[$vcount]} $var files"
         if [ $dummy -ne ${nvec[$vcount]} ]; then 
            echo "ERROR: Missing or extra $var"
            echo 24; exit 24
         fi
         mv -v ../run.*/$var* ../vectors_$it0
         ln -sf ../vectors_$it0/$var* ./
         vcount=$((vcount+1))
      done

      ex -c :"/rand_stage" +:":s/=.*/=4,/" +:wq namelist.input

      #Perform Eigen Decomp + Calculate Increment and Analysis
      mpistring="$MPICALL $DEBUGSTR -np $nproc_global $EXECUTABLE"
      #COULD MAKE THIS FASTER (MORE MEMORY) by distributing across more nodes (currently chooses first $npiens processors in $PBS_NODEFILE)

      echo "$mpistring"
      eval "$mpistring"
      if [ $SUBTIMING -eq 1 ]; then
         grep da_end_timing ../run.*/rsl.out.0000 > ../bench_time_hess-vecB.it$it0
         grep da_end_timing rsl.out.0000 > ../bench_time_finalize-riot.it$it0
      fi
   fi #rand_type == 1


#==============================================================================
   echo "Extract analysis for each outer loop, and initialize the next one..."
#==============================================================================

   if [ $it -lt $NOUTER ]; then
      if [ $it -eq $itstart ]; then cp fg fg_it0; fi
      mv -v wrfvar_output wrfvar_output_$it0
      ln -sf ./wrfvar_output_$it0 ./fg 
      ln -sf ./wrfvar_output_$it0 ./wrfinput_d01

      mkdir oldrsl_$it0
      mv -v rsl.* oldrsl_$it0
   fi
   if [ $WRF_CHEM -gt 0 ]; then
      if [ $it -eq $itstart ]; then mkdir oldhxy; fi
      mv -v *Hx_y* oldhxy/
   fi

   if [ $it -eq $itstart ]; then
      cat ../run.$ii_grad/cost_fn > ../cost_fn
      cat ../run.$ii_grad/grad_fn > ../grad_fn
   else
      # Use these lines if cost_fn file is overwritten each outer iteration
      grep [0-9] ../run.$ii_grad/cost_fn >> ../cost_fn
      grep [0-9] ../run.$ii_grad/grad_fn >> ../grad_fn
   fi

   hr0=$(($SECONDS / 3600))
   if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
   min0=$((($SECONDS / 60) % 60))
   if [ $min0 -lt 10 ]; then min0="0"$min0; fi
   sec0=$(($SECONDS % 60))
   if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
   echo "Iteration $it compute time: $hr0:$min0:$sec0"

   iteration_sec=$((iteration_sec+SECONDS))
   hr0=$(($iteration_sec / 3600))
   if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
   min0=$((($iteration_sec / 60) % 60))
   if [ $min0 -lt 10 ]; then min0="0"$min0; fi
   sec0=$(($iteration_sec % 60))
   if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
   echo "Total accum. compute time: $hr0:$min0:$sec0"

   if [ $CLEANUP -gt 0 ]; then
      cd ../
      SECONDS=0

      if [ $CLEANUP -eq 1 ]; then 
         tar -czf vectors_$it0.tar.gz vectors_$it0
         echo "Completed tar of vectors_$it0.tar.gz"
      fi

   
      rm -r vectors_$it0

      hr0=$(($SECONDS / 3600))
      if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
      min0=$((($SECONDS / 60) % 60))
      if [ $min0 -lt 10 ]; then min0="0"$min0; fi
      sec0=$(($SECONDS % 60))
      if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
      echo "Iteration $it cleanup time: $hr0:$min0:$sec0"

      cleanup_sec=$((cleanup_sec+SECONDS))
      hr0=$(($cleanup_sec / 3600))
      if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
      min0=$((($cleanup_sec / 60) % 60))
      if [ $min0 -lt 10 ]; then min0="0"$min0; fi
      sec0=$(($cleanup_sec % 60))
      if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
      echo "Total accum. cleanup time: $hr0:$min0:$sec0"

   fi
   cd $CWD
done
## END Outer Loop

SECONDS=0

# The remainder of the script should be modified to meet the needs of your application.

#Extract important output files
OUTLOC="../../$OUTDIR/$PBS_JOBNAME"
if [ $rand_type -le 6 ]; then
   mkdir ../../$OUTDIR
   mkdir $OUTLOC
   cp wrfvar_out* $OUTLOC/
   if [ $WRF_CHEM -eq 0 ]; then cp AminusB* $OUTLOC/; fi
   cp ../cost_fn $OUTLOC/
   cp ../grad_fn $OUTLOC/
   cp -r oldrsl_* $OUTLOC/
   cp rsl.* $OUTLOC/
   if [ $WRF_CHEM -gt 0 ]; then cp oldhxy/* $OUTLOC/; fi
fi

if [ $CLEANUP -gt 0 ] && [ $rand_type -ne 11 ]; then
   echo "CLEANING UP SVD $(date)"
   cd $CWD/../
   if [ $CLEANUP -eq 1 ]; then 
      mkdir oldvectors
      mv vectors_0* oldvectors/
      tar -czf oldvectors.tar.gz oldvectors
      rm -r oldvectors/

      mkdir oldrundirs
      mv run.0* oldrundirs/
      tar -czf oldrundirs.tar.gz oldrundirs
      rm -r oldrundirs/
   fi
   echo "Completed Final Cleanup"
fi

hr0=$(($SECONDS / 3600))
if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
min0=$((($SECONDS / 60) % 60))
if [ $min0 -lt 10 ]; then min0="0"$min0; fi
sec0=$(($SECONDS % 60))
if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
echo "Current cleanup time: $hr0:$min0:$sec0"

cleanup_sec=$((cleanup_sec+SECONDS))
hr0=$(($cleanup_sec / 3600))
if [ $hr0 -lt 10 ]; then hr0="0"$hr0; fi
min0=$((($cleanup_sec / 60) % 60))
if [ $min0 -lt 10 ]; then min0="0"$min0; fi
sec0=$(($cleanup_sec % 60))
if [ $sec0 -lt 10 ]; then sec0="0"$sec0; fi
echo "Total accum. cleanup time: $hr0:$min0:$sec0"

echo "FINISHED RIOT $(date)"
