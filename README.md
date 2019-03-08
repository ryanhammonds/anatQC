# filbeyQC

preprocessing:
  1. extract - bet
  2. normalize - flirt
  3. apply mask - fnirt
  
flag unusual voxels:
  1. read image to numpy array
  2. calculate standard deviation / intensity threshold of each slice
  3. flag voxels below intensity threshold that exist
  
determine pathology:
  1. count contigous unusual voxels
  2. if count is larger than sensitivity value then pathology or artifact might exist
  
