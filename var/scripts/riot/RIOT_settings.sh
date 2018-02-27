#!/bin/bash

#=======================================================================================
# INSTRUCTIONS:
#=======================================================================================
# This script defines the behavior of RIOT.  All of the below environment variables are 
# needed to run RIOT. See WRFDA_RIOT.sh for more info.
#---------------------------------------------------------------------------------------

export WRF_MET=1 #Only set > 0 for meteorology-enabled DA
export WRF_CHEM=0 #Only set > 0 for chemistry-enabled DA

export nout_RIOT=3 # number of outer iterations (overrides max_ext_its from namelist.input)
export NBLOCK=40   # number of ensemble members in first outer iteration
                   # - Similar to ntmax(1) in namelist.input, and should be larger (x2-x10)
                   #    than the CG inner iteration count to produce equivalent results.
                   # - Ensemble counts in all other outer iterations should be set in 
                   #    namelist.input as ntmax=$NBLOCK,nens[2],nens[3],nens[4], etc...
                   # - Eventually NBLOCK can be replaced with retrieval of ntmax(1) 
                   #    from namelist.input.  For now, this separate setting is a
                   #    reminder to set the number of nodes/cores.
                   # FUTURE IDEA: it may reduce wall-time of 1st ensemble member (slowest one)
                   #  to request 1 extra independent head node for process management
export rand_type=6  # 3-Block Lanczos; 6-RSVD5.6
################################################################################
#SPECIFIC TO RSVD (rand_type==6)
# - NUMNODES must be >= $((NBLOCK+1)) [where the +1 accounts for the gradient]
# - For 2 nodes per AD/TL simulation, use NUMNODES=$((2*$((NBLOCK+1)))), etc.
export prepend_rsvd_basis=0 #If ==1, prepend RSVD basis with gradient vector in 
                            # all outer iterations after the first
export ADAPT_SVD="1" # 0 (LRA), 1 (ADAPTIVE), or 2 (LRU)
export svd_p=0 # + some small value (e.g., 5) between [0,min(NBLOCK))
               # + only used when ADAPT_SVD==1
################################################################################
#SPECIFIC TO Block Lanczos (rand_type==3)
# - NUMNODES must be >= $((NBLOCK))
# - For 2 nodes per AD/TL simulation, use NUMNODES=$((2*NBLOCK)), etc.
export nin_RIOT=1  # number of Block Lanczos inner iterations
################################################################################

export GLOBAL_OPT="false" #"true" or "false"
# + GLOBAL_OPT = true
# + Global I/O is necessary to change number of processors per ensemble member between outer iterations"
# + Global cv I/O in WRFDA"
#   1 - requires additional wall-time during gather stage"
#   2 - requires additional minimum memory on head node"
#   3 - performs sequential instead of parallel I/O
# Consider these factors if da_wrfvar.exe hangs or quits unexpectedly during global cv I/O"
# Turn off global I/O with, e.g., export GLOBAL_OPT=\"false\""

export GLOBAL_OMEGA="true" #"true" or "false"
#Use GLOBAL_OMEGA="true" to ensure fully independent CV samples (OMEGA)
# + adds some sequential wall-time per outer iteration
#Use GLOBAL_OMEGA="false" if accuracy can be traded for small wall-time cost

export RIOT_RESTART=0 #If ==1, set nout_RIOT to the number of outer iterations 
                      # to complete after start file "ALT_START" 
                      #   --> posterior covariance only
                      #If ==2, set nout_RIOT to the number of outer iterations (including previous)
                      #   --> use ALT_it1 to set the alternative starting iteration
                      #   --> minimisation and posterior covariance

if [ $RIOT_RESTART -gt 0 ]; then
   #Extra files for for RIOT_RESTART>0:
   #DADIR=<location of:>
            #earlier cvt.itXX.p0000, 
            #earlier xhat.itXX.p0000, 
            #earlier wrfvar_output_XX files for linking to wrfinput_d01, fg

   PREVIOUS_RUN="SVD6_N=40_no=12"

   if [ -z "$WRFSUPER" ] && [ -z "$DADIR" ]; then #Could be defined externally
      echo ""
      echo "Need to define DADIR for RIOT_RESTART>0."
      echo "Turning off RIOT restart capability."
      echo ""
      RIOT_RESTART=0
   else
      if [ -z "$DADIR" ]; then
         #Adjust the file structure as needed
         DADIR="$WRFSUPER/DA/$PREVIOUS_RUN/run"
      fi
      export ALT_START=$DADIR"/wrfvar_output_05"

      #Three extra settings for RIOT_RESTART==2:
      export ALT_it1=9
      ii=$((ALT_it1-1))
      if [ $ii -lt 10 ]; then ii=0$ii; fi

      export ALT_CVT="$DADIR/cvt.it"$ii".p0000"
      export ALT_XHAT="$DADIR/xhat.it"$ii".p0000"
   fi
fi

export RIOT_SETTINGS_CALLED=1
