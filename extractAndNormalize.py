# flirt/bet
# 2/1/19
# Evan Kersey

# makes list of brain volumes to extract 
for i in $(ls); do
	for file in ${i}/ses-01/anat/*
	do
		if [[ $file == *.nii ]] # checks for correct file type
		then
			echo "bet ${file} ${i}/ses-01/anat/${i}_ses-01_T1w_brain.nii -R" >> runSED.txt
		fi
	done
done

# runs brain extraction with parralel processing
cat runSED.txt | /mnt/Filbey/Ryan/parallel -j 20 {}

# makes list of extracted brains to standardize
for i in $(ls); do
	for file in ${i}/ses-01/anat/*
	do
		if [[ $file == *.gz ]]
		then
			echo "/usr/local/fsl/bin/flirt -in ${file} -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain.nii.gz -out ${i}/ses-01/anat/${i}_native2standard -omat native2standard.mat -cost corratio -dof 12 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp trilinear" >> runFLIRT.txt
		fi
	done
done

# normalizes extracted brains
cat runFLIRT.txt | /mnt/Filbey/Ryan/parallel -j 20 {}
