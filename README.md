# anatQC

Flags voxels in anatomical T1w images below a defined intensity and cluster threshold to flag potentiall issues in anatomical images. Useful prior to running freesurfer.

Output:
  1. 4dClusters*.nii.gz in linear space (overlay with sub-*/native2std_lin_masked.nii.gz).
  2. 4dClusters*.tsv

Preprocessing Workflow:
  1. Brain Extraction: FSL's BET
  2. Registration: FSL's FLIRT and FNIRT
  3. Masking: FSL's applywarp standardbrain+csf mask to linear space. The linear transfomred masked image will be used in detection.
 
If these above steps produce abnormal results, the detection workflow will flag incorrectly!
 
Detection Workflow:
  1. Flag voxels below --thr standard deviations.
  2. Determine cluster via recursion.
  3. Ignore voxels adjacent to intensity of zero (outside mask).
  4. Determine cluster and remove those below --cluster threshold.
  5. Ignore clusters whose center of mass coordinates have an intensity of 0.

Processing Options:
  1. SLURM
  
To-DO:
  1. Test serial and gnu parallel options
