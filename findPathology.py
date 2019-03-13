#!/usr/bin/python3 env

# 2/18/19
# Evan Kersey
# code pulled from Venetian_Axial.py by Ryan Hammonds
# 
# to do:
#	efficiency, only parse mri once, store differences in multidiminsional array and check standard deviation from that
# 	figure out optimal sensitivity
#		find greatest difference / sd for images with no pathology and with pathology
#	

import os
import sys
import numpy as np
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=FutureWarning)
import nibabel as nib


#number of standard deviations
sensitivity = 10 


#gets nifti as input
input_img=sys.argv[1]

img = nib.load(input_img)
img_data = img.get_data()
data_obj = img.dataobj
## get dimensions
x_len = data_obj.shape[0]
y_len = data_obj.shape[1]
z_len = data_obj.shape[2]


## creates array containing the value of one standard deviation of the difference in voxel intensity for each slice
slice_dif_std = []
slice_dif_list = []
for z in range(0, z_len):
	slice_dif_list.clear() #new list for each z
	for y in range(0, y_len):
		for x in range(0, x_len):
			if img_data[x,y,z] and img_data[x+1,y,z]: # checks if both voxels exist (x direction)
				slice_dif_list.append(abs(img_data[x,y,z] - img_data[x+1,y,z])) #appends the change in intensity to array
			if img_data[x,y,z] and img_data[x,y+1,z]: # (y direction)
				slice_dif_list.append(abs(img_data[x,y,z] - img_data[x,y+1,z]))
	if len(slice_dif_list): # if the list exists for this z (the slize is not all zeros)
		slice_dif_std.append(np.std(slice_dif_list))
	else:	
		slice_dif_std.append(0) # if the slice has no voxels on it



badSlice = [0] * z_len # array of slices good or bad
for z in range(0, z_len):
	if slice_dif_std[z]: # checks if slice contains data
		for y in range(0, y_len):
			if badSlice[z]: # skips if slice has already determined to be bad
				break
			for x in range(0, x_len):
				if img_data[x,y,z] and img_data[x+1,y,z]: # if both vexels exist
					if (img_data[x,y,z] - img_data[x+1,y,z]) > (slice_dif_std[z] * sensitivity): # if slice differences are greater than the sensitivity value for slice
						badSlice[z] = True # slice is bad
						break
				if img_data[x,y,z] and img_data[x,y+1,z]:
					if (img_data[x,y,z] - img_data[x,y+1,z]) > (slice_dif_std[z] * sensitivity):
						badSlice[z] = True
						break
# find number of bad slices in image
numBadSlice = 0
for i in badSlice:
	if badSlice[i] == True:
		numBadSlice += 1

print("number of bad slices: ", numBadSlice)
