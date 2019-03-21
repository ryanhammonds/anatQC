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
Evan Kersey
Ryan Hammonds

Automated pathology detection in T1 weighted anatomical images.

Usage 

Required Arguments:
	-i,  --input	<path/to/BIDS/dir>	Parent directory containing feat directories.

	-o,  --output	<path/to/write/results>	Directory to write pathology detection results.  
							
	
	-t, --thr	<float>			Intensity standard deviations for detection.
						Ex: -t 2 will flag voxels > 2 stddev from mean voxel intensity in region.		  
    
Optional Arguments:
	-c, --cluster	<int>			Cluster size required for detection. 
						Default 10 contiguous voxels. 
	
	-s, --seg	<path/to/seg.nii>	Anatomical segmentation mask used to search signal dropout region by region. 
						Default: Harvard-Oxford Cortical and Subcortical Atlas

	-nc, --n-cpus	<int>			Number of CPUs to run processing in parallel. Ignored is --slurm is used.

	-sl, --slurm				For use on SLURM cluster. If used, REQUIRES:
						--partition
						--nodes				
						--ntasks			
						--time	<hr:min:sec>
						Note: ntasks = number sujects
						Note: May have to modify code if n subjs > (CPUs/Node x nNodes)

EOF
exit 1
}

error() {
  if [[ $REQ == ERROR ]]; then 
	echo -e "\nERROR: Missing required arguments."
  elif [[ SLURM_ERROR = 1 ]]; then
  	echo -e "\nError: Missing required slurm arguments."
  fi
  
  echo 'See --help'
  exit 1
}


while [[ $# -gt 0 ]]; do
key="$1"
case $key in
	-h|--help)
	usage
	;;
	-i|--input)
    	BIDS="$2"
    	shift; shift 
    	;;
  	-o|--output)
    	OUTPUT="$2"
    	shift; shift 
    	;;
	-t|--thr)
    	THR="$2"
    	shift; shift
		;;
	-c|--cluster)
		CLUSTER="$2"
		shift; shift
		;;
	-s|--seg)
		SEG="$2"
		shift; shift
		;;
	-nc|--ncpus)
		NCPUS="$2"
		shift; shift
		;;
	-sl|--slurm)
		SLURM=1
		shift;
		;;
	--partition)
		PART="$2"
		shift; shift
		;;
	--nodes)
		NODES="$2"
		shift; shift
		;;
	--ntasks)
		NTASKS="$2"
		shift; shift
		;;
	--time)
		TIME="$2"
		shift; shift
		;; 
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
		SLURM_ERROR = 1 && error
	elif [[ -z $NODES ]]; then
		SLURM_ERROR = 1 && error
	elif [[ -z $NTASKS ]]; then
		SLURM_ERROR = 1 && error
	elif [[ -z $TIME ]]; then
		SLURM_ERROR = 1 && error
	fi
fi

# Sets default args if none are specified
if [[ -z $CLUSTER ]] && CLUSTER=5
if [[ -z $SEG ]] && SEG=$FSLDIR/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr50-2mm.nii.gz


# Makes subj folders in output dir
for i in $BIDS/sub*; do 
	subj=$(echo $i | sed "s#.*/##")
	mkdir $OUTPUT/$subj
done

mkdir $OUTPUT/scripts

# Checks for multi or single session BIDS and stores anatomical paths
anatPATHS=()
firstSubj=$(ls -1 $BIDS | head -n 1)
if [[ -d $BIDS/$firstSubj/ses-01 ]]; then
	anatPATHS+=(`find $BIDS/sub*/ses-*/anat/*_T1w.nii`)
else
	anatPATHS+=(`find $BIDS/sub*/anat/*_T1w.nii`)
fi

## Runs BET, FLIRT, FNIRT ##
# Serial
if [[ -z $SLURM && -z $NCPUS ]]; then
	for image in ${anatPATHS[@]}; do
		subj=$(echo $image | sed "s/.*sub/sub/" | sed "s/_.*//")
		$FSLDIR/bin/bet $image -R
		$FSLDIR/bin/flirt -in "$image"_brain -ref $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz -out $OUTPUT/$subj/native2standard -omat $OUTPUT/$subj/native2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear
		$FSLDIR/bin/fnirt --iout=$OUTPUT/$subj/native2standard_nonlin --in=/$OUTPUT/$subj/*_brain.nii.gz --aff=$OUTPUT/$subj/native2standard.mat --cout=$OUTPUT/$subj/native2standard_warp --iout=$OUTPUT/$subj/native2standard_nonlin --jout=/$OUTPUT/$subj/native2native_jac --config=T1_2_MNI152_2mm --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz --refmask=$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil --warpres=10,10,10	
	done
# SLURM
elif [[ -n $SLURM ]]; then
	echo '#!/bin/bash'										>  $OUTPUT/scripts/slurmBET.sh
	echo \#SBATCH --partition="$PART"		 				>> $OUTPUT/scripts/slurmBET.sh
	echo \#SBATCH --nodes="$NODES" 							>> $OUTPUT/scripts/slurmBET.sh
	echo \#SBATCH --ntasks=$NTASKS							>> $OUTPUT/scripts/slurmBET.sh
	echo \#SBATCH --cpus-per-task=1			 				>> $OUTPUT/scripts/slurmBET.sh
	echo \#SBATCH --time="$TIME"							>> $OUTPUT/scripts/slurmBET.sh
	echo -e "\n"											>> $OUTPUT/scripts/slurmBET.sh
	for image in ${anatPATHS[@]}; do
		echo "srun -N 1 -n 1 $FSLDIR/bin/bet $image -R" 	>> $OUTPUT/scripts/slurmBET.sh 
	done
	echo 'wait'												>> $OUTPUT/scripts/slurmBET.sh	
	sbatch $OUTPUT/scripts/slurmBET.sh
	#FLIRT
	cat $OUTPUT/scripts/slurmBET.sh | head -n 7				> $OUTPUT/scripts/slurmFLIRT.sh
	for image in ${anatPATHS[@]}; do
		subj=$(echo $image | sed "s/.*sub/sub/" | sed "s/_.*//")
		echo "srun -N 1 -n 1 $FSLDIR/bin/flirt -in "$image"_brain -ref $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz -out $OUTPUT/$subj/native2standard -omat $OUTPUT/$subj/native2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear" >> $OUTPUT/scripts/slurmFLIRT.sh 
		echo "wait"											>> $OUTPUT/scripts/slurmFLIRT.sh
	done
	sbatch $OUTPUT/scripts/slurmFLIRT.sh
	#FNIRT
	cat $OUTPUT/scripts/slurmBET.sh | head -n 7				> $OUTPUT/scripts/slurmFNIRT.sh
	for image in ${anatPATHS[@]}; do
		subj=$(echo $image | sed "s/.*sub/sub/" | sed "s/_.*//")
		echo "srun -N 1 -n 1 $FSLDIR/bin/fnirt --iout=$OUTPUT/$subj/native2standard_nonlin --in=/$OUTPUT/$subj/*_brain.nii.gz --aff=$OUTPUT/$subj/native2standard.mat --cout=$OUTPUT/$subj/native2standard_warp --iout=$OUTPUT/$subj/native2standard_nonlin --jout=/$OUTPUT/$subj/native2native_jac --config=T1_2_MNI152_2mm --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz --refmask=$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil --warpres=10,10,10" >> $OUTPUT/scripts/slurmFNIRT.sh
		echo "wait"											>> $OUTPUT/scripts/slurmFNIRT.sh
	done
	sbatch $OUTPUT/scripts/slurmFNIRT.sh
# GNU Parallel
elif [[ -n $NCPUS && -z $SLURM ]]; then
	for image in ${anatPATHS[@]}; do
		subj=$(echo $image | sed "s/.*sub/sub/" | sed "s/_.*//")
		echo "$FSLDIR/bin/bet $image -R" >> $OUTPUT/runBet 
		echo "$FSLDIR/bin/flirt -in "$image"_brain -ref $FSLDIR/data/standard/MNI152_T1_1mm_brain.nii.gz -out $OUTPUT/$subj/native2standard -omat $OUTPUT/$subj/native2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear" >> $OUTPUT/scripts/runFLIRT
		echo "fnirt --iout=$OUTPUT/$subj/native2standard_nonlin --in=/$OUTPUT/$subj/*_brain.nii.gz --aff=$OUTPUT/$subj/native2standard.mat --cout=$OUTPUT/$subj/native2standard_warp --iout=$OUTPUT/$subj/native2standard_nonlin --jout=/$OUTPUT/$subj/native2native_jac --config=T1_2_MNI152_2mm --ref=$FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz --refmask=$FSLDIR/data/standard/MNI152_T1_2mm_brain_mask_dil --warpres=10,10,10" >> $OUTPUT/scripts/runFNIRT	
	done
	baseDIR=$(echo $0 | sed "s#/.*#/#")
	cat $OUTPUT/scripts/runBET | $baseDIR/parallel -j $NCPUS
	cat $OUTPUT/scripts/runFLIRT | $baseDIR/parallel -j $NCPUS
	cat $OUTPUT/scripts/runFNIRT | $baseDIR/parallel -j $NCPUS
fi
wait

## Call python script here ##
# 1) Will need to adjust python script to pull images from $OUTPUT/$subj/native2standard_nonlin
# 2) Arguments set in this script can carry over to findPathology.py
#    Ex: python3.6 findPathology.py $THR
#  		 Then in findPathology.py $THR become sys.argv[1]
#		 Note: global python3.6 will not have packages you need (nibabel, etc). You will need to change it to environment python3

