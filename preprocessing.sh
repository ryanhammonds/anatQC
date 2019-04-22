#!/bin/bash

## TO-DO ##
# 1) Add check that verify BIDS dir is properly structured and named.
# 2) Call pathology.py in this script - look at the bootom of this script for what's next to do.
# 3) This script assumes pathology.py parses anatomical regions one by one to look for low intensity clusters in comparison
#	 to the mean of the region the cluster is in. This will be an issue if there is a hole that spans an entire region and
#	 the mean of that region becomes very low. Maybe the mean value to compare voxels too (within single region) should come
#  	 from a mean of region across all participants.	There are other ways to solve this though. Do what you think is most efficient.
#  	 I haven't made it to your findPathology.py script yet so I'm not sure what approach you've decided on.

usage() {
  cat << EOF >&2

anatQC

Authors:
Ryan Hammonds
Evan Kersey

Automated cluster signal dropout detection in T1 weighted anatomical images.

Usage

Required Arguments:
	-i,  --input	<path/to/BIDS/dir>	Parent directory containing BIDS directories.

	-o,  --output	<path/to/write/results>	Directory to write pathology detection results.

Optional Arguments:
	-t, --thr	<float>			Intensity standard deviations for detection.
						Defaults is 2 standard deviations.

						Note: if the --thr produces a negative voxel intensity, --thr will be
						reduce in increments of 0.1 until --thr 1 or until a positive intensity
						is reached.

	-c, --cluster	<int>			Cluster size required for detection.
						Default 50 contiguous voxels.

	-nc, --n-cpus	<int>			Number of CPUs to run processing using GNU parallel.
						Ignored if --slurm is used. Use --ntasks-per-node instead.
						Default: Serial processing on one CPU.

SLURM Arguments:
	--slurm			No argmuent. 	Indicates argments below are required.
	--partition		<debug>		Name of paritition.
	--nodes			<int>		Total number of nodes.
	--ntasks-per-node	<int>		Max number of subjects/images to run on each nodes.
	--time			<hr:min:sec>	Time limit.

EOF
exit 1
}

error() {
  if [[ $REQ == ERROR ]]; then
    echo -e "\nMissing required arguments."
  elif [[ SLURM_ERROR = 1 ]]; then
    echo -e "\nMissing required slurm arguments."
  fi

  if [[ $NON_ARG == ERROR ]]; then
    echo -e "\nUnrecognized argument: $NON_ARG_VAL"
  fi

  echo -e 'See --help\n'
  exit 1
}


while [[ $# -gt 0 ]]; do
key="$1"
case $key in
	-h|--help)
	usage
	;;
	-i|--input)
    	BIDS=$2
    	shift; shift
    	;;
  	-o|--output)
    	OUTPUT=$2
    	shift; shift
    	;;
	-t|--thr)
    	THR=$2
    	shift; shift
		;;
	-c|--cluster)
		CLUSTER=$2
		shift; shift
		;;
	-nc|--ncpus)
		NCPUS=$2
		shift; shift
		;;
	-sl|--slurm)
		SLURM=1
		shift
		;;
	--partition)
		PART=$2
		shift; shift
		;;
	--nodes)
		NODES=$2
		shift; shift
		;;
	--ntasks-per-node)
		NTASKS=$2
		shift; shift
		;;
	--time)
		TIME=$2
		shift; shift
		;;
	-*|--*)
		NON_ARG=ERROR
		NON_ARG_VAL=$2
		error
		;;
	*)
		break
esac
done


# Checks if required arg are missing
[[ -z $BIDS ]] 	&& REQ=ERROR
[[ -z $OUTPUT ]] 	&& REQ=ERROR
[[ -z $THR ]] 		&& REQ=ERROR
[[ $REQ == ERROR ]]	&& error

# Checks SLURM and its arguments
if [[ -n $SLURM ]]; then
	if [[ -z $PART ]]; then
		SLURM_ERROR=1 && error
	elif [[ -z $NODES ]]; then
		SLURM_ERROR=1 && error
	elif [[ -z $NTASKS ]]; then
		SLURM_ERROR=1 && error
	elif [[ -z $TIME ]]; then
		SLURM_ERROR=1 && error
	fi
fi

baseDir=$(echo $0 | sed "s/\/preprocessing.sh//")

# Sets default args if none are specified
[[ -z $CLUSTER ]] && CLUSTER=50
[[ -z $THR ]] && THR=2

# Remove trailing slash
OUTPUT=${OUTPUT%/}
BIDS=${BIDS%/}

[[ ! -d $OUTPUT ]] && mkdir $OUTPUT
[[ ! -d $OUTPUT/scripts ]] && mkdir $OUTPUT/scripts

# Checks for multi or single session BIDS and stores anatomical paths
anatPATHS=()
firstSubj=$(ls -1 $BIDS | grep -e "^sub.*" | head -n 1)
if [[ -d $BIDS/$firstSubj/ses-01 ]]; then
	anatPATHS+=(`find $BIDS/sub*/ses-*/anat/*_T1w.nii`)
else
	anatPATHS+=(`find $BIDS/sub*/anat/*_T1w.nii`)
fi

# Make subject folders for those who have T1w images.
for i in ${anatPATHS[@]}; do
	subj=$(echo $i | sed "s#.*/##" | sed "s/_.*//")
	[[ ! -d $OUTPUT/$subj ]] && mkdir $OUTPUT/$subj
	# Remove QC.py cluster outputs if they exist
	[[ -f $OUTPUT/$subj/4dClusters.nii.gz ]] && rm $OUTPUT/$subj/4dClusters.nii.gz
	[[ -f $OUTPUT/$subj/4dClusters.tsv ]] && rm $OUTPUT/$subj/4dClusters.tsv
	[[ -f $OUTPUT/$subj/4dClusters.tsv ]] && rm $OUTPUT/$subj/voxel_thr_used.txt
done

# Make text file for SLURM access
[[ -f $OUTPUT/scripts/$image/image_paths.txt ]] && rm $OUTPUT/scripts/$image/image_paths.txt
[[ -f $OUTPUT/flagged_ids.txt ]] && rm $OUTPUT/flagged_ids.txt
for image in ${anatPATHS[@]}; do
	echo $image	>> $OUTPUT/scripts/image_paths.txt
done

## Runs BET, FLIRT, FNIRT ##
# Serial
if [[ -z $SLURM && -z $NCPUS ]]; then
	for image in ${anatPATHS[@]}; do
		subj=$(echo $image | sed "s/.*sub/sub/" | sed "s/_.*//")
		T1w=$(echo $image | sed "s/.*sub/sub/" | sed "s/\.nii//g")
		$FSLDIR/bin/bet $image $OUTPUT/$subj/"$T1w"_brain -R
		$FSLDIR/bin/flirt -in $OUTPUT/$subj/"$T1w"_brain -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -out $OUTPUT/$subj/native2std -omat $OUTPUT/$subj/native2std.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear &

		# Transform standard space masks into linear MNI space.
		$FSLDIR/bin/convert_xfm -omat $OUTPUT/$subj/native2std_inv.mat -inverse $OUTPUT/$subj/native2std.mat
		$FSLDIR/bin/flirt -in $baseDir/masks/mask.nii.gz -ref $T1w -out $OUTPUT/$subj/mask_native -applyxfm -init $OUTPUT/$subj/native2std_inv.mat
		$FSLDIR/bin/fslmaths $OUTPUT/$subj/mask_native.nii.gz -bin $OUTPUT/$subj/mask_native.nii.gz
		$FSLDIR/bin/fslmaths $OUTPUT/$subj/native2std.nii.gz -mas $OUTPUT/$subj/mask_native $OUTPUT/$subj/native2std_masked
		done
# SLURM
elif [[ -n $SLURM ]]; then

	# Setup SLURM logs
	if [[ -d $OUTPUT/scripts/slurmLogs ]]; then
		rm $OUTPUT/scripts/squegs/*
	else
		mkdir $OUTPUT/scripts/slurmLogs
	fi

  [[ ! -d $OUTPUT/scripts/wf_logs ]] && mkdir $OUTPUT/scripts/wf_logs

	SUBJS=${#anatPATHS[@]}
	CPUS=$(($NTASKS*$NODES))

  echo '#!/bin/bash'														                           >  $OUTPUT/scripts/slurm_preproc_wf.sh
  echo \#SBATCH --partition=$PART		 								                       >> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo \#SBATCH --cpus-per-task=1			 							                     	 >> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo \#SBATCH --time=$TIME											                      	 >> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo \#SBATCH --array=0-"$SUBJS"%$CPUS								                 	 >> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo \#SBATCH --output=$OUTPUT/scripts/slurmLogs/preproc_wf_%A_%4a.log   >> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo -e "\n"															                               >> $OUTPUT/scripts/slurm_preproc_wf.sh
  # Create array to reference for anatomical inputs, output directory, and subject
	echo 'anatArray=(); for image in $(cat $OUTPUT/scripts/image_paths.txt); do anatArray+=($image); done'	>> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo 'outArray=(); idArray=(); for subjDir in $(ls -d1 $OUTPUT/sub-*); do out=$(echo $subjDir); outArray+=($out); id=$(echo $subjDir | sed "s#.*/##"); idArray+=($id); done' >> $OUTPUT/scripts/slurm_preproc_wf.sh

	# Create workflow, then export it to be accessed from sbatch
  preproc_wf() {
		# Check if steps have been executed
    if [[ -f $OUTPUT/scripts/wf_logs/"$3".log ]]; then
      reorientCheck=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'Reorient Complete')
      betCheck=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'BET Complete')
      flirtCheck=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'FLIRT Complete')
      fnirtCheck=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'FNIRT Complete')
      apply1Check=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'Apply1 Complete')
      invCheck=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'Inv Complete')
      apply2Check=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'Apply2 Complete')
      apply3Check=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'Apply3 Complete')
      math1Check=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'Math1 Complete')
      math2Check=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'Math2 Complete')
      #QC_Check=$(cat $OUTPUT/scripts/wf_logs/"$3".log | grep 'QC_Check Complete')
    fi

    # Reorientation
    if [[ -z $reorientCheck ]]; then
        $FSLDIR/bin/fslreorient2std $1 $2/"$3"_T1w
        wait
        echo 'Reorient Complete' >> $OUTPUT/scripts/wf_logs/"$3".log
    fi
    # Brain extraction
    if [[ -z $betCheck ]]; then
      $FSLDIR/bin/bet $2/"$3"_T1w $2/"$3"_T1w_brain -R
		  wait
      echo 'BET Complete' >> $OUTPUT/scripts/wf_logs/"$3".log
    fi
		# Linear registration on brain extracted T1w
    if [[ -z $flirtCheck ]]; then
  		$FSLDIR/bin/flirt -in $2/"$3"_T1w_brain \
      -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
      -out $2/native2std_lin -omat $2/native2std_lin.mat \
      -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear
      wait
      echo 'FLIRT Complete' >> $OUTPUT/scripts/wf_logs/"$3".log
    fi
		# Non-linear registration on whole T1w.
    if [[ -z $fnirtCheck ]]; then
		    $FSLDIR/bin/fnirt --iout=$2/native2std_nonlin_head \
        --in=$2/"$3"_T1w \
        --aff=$2/native2std_lin.mat \
        --cout=$2/native2std_nonlin_warp \
        --iout=$2/native2std_nonlin \
        --jout=$2/native2native_nonlin_jac \
        --config=T1_2_MNI152_2mm \
        --ref=$FSLDIR/data/standard/MNI152_T1_2mm \
        --refmask=$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask \
        --warpres=10,10,10
        wait
        echo 'FNIRT Complete' >> $OUTPUT/scripts/wf_logs/"$3".log
    fi
		# Apply above transformation to brain extracted image.
    if [[ -z $applyCheck ]]; then
		    $FSLDIR/bin/applywarp -i $2/"$3"_T1w_brain \
        -r $FSLDIR/data/standard/MNI152_T1_2mm_brain \
        -o $2/native2std_nonlin_brain \
        -w $2/native2std_nonlin_warp
        wait
        echo 'Apply1 Complete' >> $OUTPUT/scripts/wf_logs/"$3".log
    fi
		# Inverse the warp field for standard to native space transformations
    if [[ -z $invCheck ]]; then
      $FSLDIR/bin/invwarp \
      --ref=$2/"$3"_T1w_brain \
      --warp=$2/native2std_nonlin_warp \
      --out=$2/native2std_nonlin_warp_inv
  		wait
      echo 'Inv Complete' >> $OUTPUT/scripts/wf_logs/"$3".log
    fi
		# Warp brain+CSF mask into native space
		if [[ -z $apply2Check ]]; then
      $FSLDIR/bin/applywarp \
      -i $baseDir/masks/mask_2mm.nii.gz \
      -o $2/mask_native -r $2/"$3"_T1w_brain \
      --warp=$2/native2std_nonlin_warp_inv
  		wait
      echo 'Apply2 Complete' >> $OUTPUT/scripts/wf_logs/"$3".log
    fi
		# Warp native space mask into linear space
		if [[ -z $apply3Check ]]; then
      $FSLDIR/bin/applywarp \
      -i $2/mask_native \
      -r $2/native2std_lin.nii.gz \
      --premat=$2/native2std_lin.mat \
      -o $2/mask_lin
  		wait
      echo 'Apply3 Complete' >> $OUTPUT/scripts/wf_logs/"$3".log
    fi
		# Linear image masking
    if [[ -z $math1Check ]]; then
		    $FSLDIR/bin/fslmaths $2/mask_lin -bin $2/mask_lin
        wait
        echo 'Math1 Complete' >> $OUTPUT/scripts/wf_logs/"$3".log
    fi
    if [[ -z $math2Check ]]; then
		    $FSLDIR/bin/fslmaths $2/native2std_lin.nii.gz -mas $2/mask_lin $2/native2std_lin_masked
        wait
        echo 'Math2 Complete' >> $OUTPUT/scripts/wf_logs/"$3".log
    fi
		# Run py script on linear image
    $baseDir/QC.py $2/native2std_lin_masked.nii.gz $THR $CLUSTER
    wait
        #echo 'QC_Check Complete' >> $OUTPUT/scripts/wf_logs/"$3".log

	}
	# Export for SLURM access.
	export -f preproc_wf
	export OUTPUT
	export CLUSTER
	export THR
	export baseDir

  echo '' >> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo 'preproc_wf ${anatArray[SLURM_ARRAY_TASK_ID]} ${outArray[SLURM_ARRAY_TASK_ID]} ${idArray[SLURM_ARRAY_TASK_ID]}'	>> $OUTPUT/scripts/slurm_preproc_wf.sh
	sed -i "0,/\/\//s//\//" $OUTPUT/scripts/slurm_preproc_wf.sh
	sbatch $OUTPUT/scripts/slurm_preproc_wf.sh

# GNU Parallel
elif [[ -n $NCPUS && -z $SLURM ]]; then
	for image in ${anatPATHS[@]}; do
		subj=$(echo $image | sed "s/.*sub/sub/" | sed "s/_.*//")
		T1w=$(echo $image | sed "s/.*sub/sub/" | sed "s/\.nii//g")
		echo "$FSLDIR/bin/bet $image $OUTPUT/$subj/"$T1w"_brain -R &" 	>> $OUTPUT/scripts/gnuBET
		echo "$FSLDIR/bin/flirt -in $OUTPUT/$subj/"$T1w"_brain -ref $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz -out $OUTPUT/$subj/native2std -omat $OUTPUT/$subj/native2std.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear &" >> $OUTPUT/scripts/gnuFLIRT
		echo "$FSLDIR/bin/convert_xfm -omat $OUTPUT/$subj/native2std_inv.mat -inverse $OUTPUT/$subj/native2std.mat &"	>> $OUTPUT/scripts/gnuTransXFM
		echo "$FSLDIR/bin/flirt -in masks/mask.nii.gz -ref $T1w -out $OUTPUT/$subj/mask_native -applyxfm -init $OUTPUT/$subj/native2std_inv.mat &"	>> $OUTPUT/scripts/gnuTransFLIRT
		echo "$FSLDIR/bin/fslmaths $OUTPUT/$subj/mask_native.nii.gz -bin $OUTPUT/$subj/mask_native.nii.gz &"	>> $OUTPUT/scripts/gnuTransMATHS
		echo "$FSLDIR/bin/fslmaths $OUTPUT/$subj/native2std.nii.gz -mas $OUTPUT/$subj/mask_native $OUTPUT/$subj/native2std_masked &" >> $OUTPUT/scripts/gnuTransMATHS2
	done

	cat $OUTPUT/scripts/gnuBET | ./parallel -j $NCPUS
	cat $OUTPUT/scripts/gnuFLIRT | ./parallel -j $NCPUS
	cat $OUTPUT/scripts/gnuTransXFM | ./parallel -j $NCPUS
	cat $OUTPUT/scripts/gnuTransFLIRT | ./parallel -j $NCPUS
	cat $OUTPUT/scripts/gnuTransMATHS | ./parallel -j $NCPUS
	cat $OUTPUT/scripts/gnuTransMATHS2 | ./parallel -j $NCPUS
fi

for clust_img in $(ls -1 $OUTPUT/sub*/4d*.nii.gz); do 
  id=$(echo $clust_img | head -n 1 | sed "s/.*sub/sub/" | sed "s#/.*##")
  echo $id >> $OUTPUT/flagged_ids.txt
done

## Need to figure out how to tie in req packages in repo with needing to source env
#source ~/env_3.6/bin/activate

#IMG=$OUTPUT/$subj/native2std_nonlin.nii.gz

#./findPathology.py $IMG $THR $CLUSTER

## Call python script here ##
# 1) Will need to adjust python script to pull images from $OUTPUT/$subj/native2std_nonlin
# 2) Arguments set in this script can carry over to findPathology.py
#    Ex: python3.6 findPathology.py $THR
#  		 Then in findPathology.py $THR become sys.argv[1]
#		 Note: global python3.6 will not have packages you need (nibabel, etc). You will need to change it to environment python3
