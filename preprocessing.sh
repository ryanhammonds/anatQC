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

Authors:
Ryan Hammonds
Evan Kersey

Automated pathology detection in T1 weighted anatomical images.

Usage 

Required Arguments:
	-i,  --input	<path/to/BIDS/dir>	Parent directory containing feat directories.

	-o,  --output	<path/to/write/results>	Directory to write pathology detection results.  
							
	
	-t, --thr	<float>			Intensity standard deviations for detection.
						Ex: -t 2 will flag voxels > 2 stddev from mean voxel intensity in region.		  
    
Optional Arguments:
	-c, --cluster	<int>			Cluster size required for detection. 
						Default 50 contiguous voxels. 

	-nc, --n-cpus	<int>			Number of CPUs to run processing in parallel. Ignored if --slurm is used. Use --ntasks instead.

	-sl, --slurm				For use on SLURM cluster. If used, REQUIRES:
						--partition
						--nodes				
						--ntasks-per-node	<int>		Max number of subjects/images to run on each nodes.			
						--time	<hr:min:sec>
						--subjs

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

# Sets default args if none are specified
[[ -z $CLUSTER ]] && CLUSTER=50

# Makes subj folders in output dir
for i in $BIDS/sub*; do 
	subj=$(echo $i | sed "s#.*/##")
	[[ ! -d $OUTPUT/$subj ]] && mkdir $OUTPUT/$subj
done

[[ ! -d $OUTPUT/scripts ]] && mkdir $OUTPUT/scripts

# Checks for multi or single session BIDS and stores anatomical paths
anatPATHS=()
firstSubj=$(ls -1 $BIDS | grep -e "^sub.*" | head -n 1)
if [[ -d $BIDS/$firstSubj/ses-01 ]]; then
	anatPATHS+=(`find $BIDS/sub*/ses-*/anat/*_T1w.nii`)
else
	anatPATHS+=(`find $BIDS/sub*/anat/*_T1w.nii`)
fi

# Bash doesn't support array export, so save to text file to read into sbatch script.
[[ -f $OUTPUT/scripts/$image/image_paths.txt ]] && rm $OUTPUT/scripts/$image/image_paths.txt
for image in ${anatPATHS[@]}; do
	echo $image	>> $OUTPUT/scripts/image_paths.txt
done

# Need this for masks paths
baseDir=$(echo $0 | sed "s/\/preprocessing.sh//")
export baseDir

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
	[[ ! -d $OUTPUT/scripts/slurmLogs ]] && mkdir $OUTPUT/scripts/slurmLogs
	SUBJS=${#anatPATHS[@]}
	CPUS=$(($NTASKS*$NODES))

	echo '#!/bin/bash'														>  $OUTPUT/scripts/slurm_preproc_wf.sh
	echo \#SBATCH --partition=$PART		 									>> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo \#SBATCH --cpus-per-task=1			 								>> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo \#SBATCH --time=$TIME												>> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo \#SBATCH --array=0-"$SUBJS"%$CPUS									>> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo \#SBATCH --output=$OUTPUT/scripts/slurmLogs/preproc_wf_%A_%4a.log	>> $OUTPUT/scripts/slurm_preproc_wf.sh		
	echo -e "\n"															>> $OUTPUT/scripts/slurm_preproc_wf.sh
	# Create array to reference for anatomical inputs, output directory, and subject
	echo 'anatArray=(); for image in $(cat $OUTPUT/scripts/image_paths.txt); do anatArray+=($image); done'	>> $OUTPUT/scripts/slurm_preproc_wf.sh
	echo 'outArray=(); idArray=(); for subjDir in $(ls -d1 $OUTPUT/sub-*); do out=$(echo $subjDir | sed "s/\/\//\//"); outArray+=($out); id=$(echo $subjDir | sed "s/\/\//\//" | sed "s/.*\///"); idArray+=($id); done' >> $OUTPUT/scripts/slurm_preproc_wf.sh

	# Create workflow, then export it to be accessed from sbatch
	preproc_wf() {
		# Brain extraction
		$FSLDIR/bin/bet $1 $2/"$3"_T1w_brain -R
		wait
		# Linear registration on brain extracted T1w
		$FSLDIR/bin/flirt -in $2/"$3"_T1w_brain -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz -out $2/native2std_lin -omat $2/native2std_lin.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear
		wait
		# Non-linear registration on whole T1w.		
		$FSLDIR/bin/fnirt --iout=$2/native2std_nonlin_head --in=$1 --aff=$2/native2std_lin.mat --cout=$2/native2std_nonlin_warp --iout=$2/native2std_nonlin --jout=$2/native2native_nonlin_jac --config=T1_2_MNI152_2mm --ref=$FSLDIR/data/standard/MNI152_T1_2mm --refmask=$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask --warpres=10,10,10
		wait
		# Apply above transformation to brain extracted image.	
		$FSLDIR/bin/applywarp -i $1 -r $FSLDIR/data/standard/MNI152_T1_2mm_brain -o $2/native2std_nonlin_brain -w $2/native2std_nonlin_warp
		wait		
		# Inverse the warp field for standard to native space transformations		
		$FSLDIR/bin/invwarp --ref=$2/"$3"_T1w_brain --warp=$2/native2std_nonlin_warp --out=$2/native2std_nonlin_warp_inv			
		wait		
		# Warp brain+CSF mask into native space
		$FSLDIR/bin/applywarp -i $baseDir/masks/mask_2mm.nii.gz -o $2/mask_native -r $2/"$3"_T1w_brain --warp=$2/native2std_nonlin_warp_inv		
		wait		
		# Warp native space mask into linear space		
		$FSLDIR/bin/applywarp -i $2/mask_native -r $2/native2std_lin.nii.gz --premat=$2/native2std.mat -o $2/mask_lin
		wait
		# Linear image masking
		$FSLDIR/bin/fslmaths $2/mask_lin -bin $2/mask_lin
		wait
		$FSLDIR/bin/fslmaths $2/native2std_lin.nii.gz -mas $2/mask_lin $2/native2std_lin_masked
		wait
		
		# Run py script on linear image
		
	}
	export -f preproc_wf
	export OUTPUT
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

