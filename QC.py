#!/usr/bin/env python

import sys
import re
import time

import numpy as np
import nibabel as nib
from scipy.ndimage import center_of_mass
from scipy.ndimage import label

# gets nifti as input
input_img=sys.argv[1]
voxel_std = float(sys.argv[2])
clust_size = int(sys.argv[3])

#load nifti to object
img = nib.load(input_img)
img_data = img.get_data()
data_obj = img.dataobj

# get dimensions
x_len = data_obj.shape[0]
y_len = data_obj.shape[1]
z_len = data_obj.shape[2]

# creates a list of all nonzero intensities, used to calculate std/mean
voxelValsList = []
for z in range(0, z_len):
	for y in range(0, y_len):
		for x in range(0, x_len):
			if img_data[x,y,z] > 0:
				voxelValsList.append(img_data[x,y,z])

meanIntensity = np.mean(voxelValsList)
stdIntensity = np.std(voxelValsList)
threshold = meanIntensity - (voxel_std * stdIntensity)

# list of all coordinates below threshold
flagList = []
for z in range(0, z_len):
	for y in range(0, y_len):
		for x in range(0, x_len):
			if img_data[x,y,z]:
				if img_data[x,y,z] < threshold:
					flagList.append([x,y,z])

# recursive function, finds clusters in 3D 
sys.setrecursionlimit(len(flagList))

def checkContig(pt, flagList, clustList):
	x = pt[0]
	y = pt[1]
	z = pt[2]

	if pt in flagList:
		clustList.append(pt)
		flagList.remove(pt)

	# All voxels surrounding a given voxels (ie rubic's cube)
	points = [[x,y,z-1], [x+1,y,z-1], [x+1,y+1,z-1], [x+1,y-1,z-1], [x-1,y,z-1], [x-1,y+1,z-1], [x-1,y-1,z-1], [x,y+1,z-1], [x,y-1,z-1], [x,y,z+1], [x+1,y,z+1], [x+1,y+1,z+1], [x+1,y-1,z+1], [x-1,y,z+1], [x-1,y+1,z+1], [x-1,y-1,z+1], [x,y+1,z+1], [x,y-1,z+1], [x,y-1,z], [x+1,y-1,z], [x+1,y-1,z+1], [x+1,y-1,z-1], [x-1,y-1,z], [x-1,y-1,z+1], [x-1,y-1,z-1], [x,y-1,z+1], [x,y-1,z-1], [x,y+1,z], [x+1,y+1,z], [x+1,y+1,z+1], [x+1,y+1,z-1], [x-1,y+1,z], [x-1,y+1,z+1], [x-1,y+1,z-1], [x,y+1,z+1], [x,y+1,z-1], [x-1,y,z], [x-1,y+1,z], [x-1,x+1,z+1], [x-1,y+1,z-1], [x-1,y-1,z], [x-1,y-1,z+1], [x-1,y-1,z-1], [x-1,y,z+1], [x-1,y,z-1], [x+1,y,z], [x+1,y+1,z], [x+1,y+1,z+1], [x+1,y+1,z-1], [x+1,y-1,z], [x+1,y-1,z+1], [x+1,y-1,z-1], [x+1,y,z+1], [x+1,y,z-1]]
	
	# Remove coords adjacent to at least 2 voxels w/intensity of 0. 
	boundary_zeros = []
	for coord in points:
		if img_data[coord[0]][coord[1]][coord[2]] == 0:
			boundary_zeros.append(0)
		if len(boundary_zeros) > 2:
			clustList.remove(pt)
			break
	
	for coord in points:
		if coord in flagList and coord not in clustList:
			clustList.append(coord)	
			flagList.remove(coord)
			checkContig(coord, flagList, clustList)
	return clustList


def checkBoundary(cluster):
	# Checks if cluster is a ring around brain boundary
	boundary_thresh = len(cluster) * 2
	boundary_zeros = []
	for voxel in cluster:
		x = voxel[0]
		y = voxel[1]
		z = voxel[2]	
		points = [[x,y,z-1], [x+1,y,z-1], [x+1,y+1,z-1], [x+1,y-1,z-1], [x-1,y,z-1], [x-1,y+1,z-1], [x-1,y-1,z-1], [x,y+1,z-1], [x,y-1,z-1], [x,y,z+1], [x+1,y,z+1], [x+1,y+1,z+1], [x+1,y-1,z+1], [x-1,y,z+1], [x-1,y+1,z+1], [x-1,y-1,z+1], [x,y+1,z+1], [x,y-1,z+1], [x,y-1,z], [x+1,y-1,z], [x+1,y-1,z+1], [x+1,y-1,z-1], [x-1,y-1,z], [x-1,y-1,z+1], [x-1,y-1,z-1], [x,y-1,z+1], [x,y-1,z-1], [x,y+1,z], [x+1,y+1,z], [x+1,y+1,z+1], [x+1,y+1,z-1], [x-1,y+1,z], [x-1,y+1,z+1], [x-1,y+1,z-1], [x,y+1,z+1], [x,y+1,z-1], [x-1,y,z], [x-1,y+1,z], [x-1,x+1,z+1], [x-1,y+1,z-1], [x-1,y-1,z], [x-1,y-1,z+1], [x-1,y-1,z-1], [x-1,y,z+1], [x-1,y,z-1], [x+1,y,z], [x+1,y+1,z], [x+1,y+1,z+1], [x+1,y+1,z-1], [x+1,y-1,z], [x+1,y-1,z+1], [x+1,y-1,z-1], [x+1,y,z+1], [x+1,y,z-1]]	
		for adjacent_voxel in points:
			if img_data[adjacent_voxel[0]][adjacent_voxel[1]][adjacent_voxel[2]] == 0:
				boundary_zeros.append(0)
				if len(boundary_zeros) > boundary_thresh:
					return False
					break
				else:
					return True		


# flag clusters
clusters = []
for xyz in flagList:
	singleCluster = []	
	singleCluster = checkContig(xyz, flagList, singleCluster)
	clusters.append(singleCluster)

# remove boundary clusters
cleaned_clusters = []
for cluster in clusters:
	if len(cluster) > clust_size:
		#if checkBoundary(cluster):
		cleaned_clusters.append(cluster)


if len(cleaned_clusters) > 0:
	# create mask of clusters and output sizes
	for idx,cluster in enumerate(cleaned_clusters):
		mask = np.zeros(np.shape(img_data))
		for voxel in cluster:
				mask[voxel[0], voxel[1], voxel[2]] = 1
		if idx == 0:	
			mask4D = mask[..., np.newaxis]			
		else:
			mask[voxel[0], voxel[1], voxel[2]] = 1
			mask = mask[..., np.newaxis]
			mask4D = np.concatenate((mask4D,mask), axis=3)

centers = []		
for cluster_idx in range(0, len(cleaned_clusters)):
	center = center_of_mass(mask4D[:,:,:,cluster_idx])
	list(center)	
	center = [int(round(dist)) for dist in center]
	centers.append(center)
	
'''
# calculate center of mass for each cluster
if len(cleaned_clusters) > 1:
	clust_len = [val for val in range(1,len(cleaned_clusters)+1)]
	lbl = label(mask)[0]
	centers = center_of_mass(mask, lbl, clust_len)
	for centerIdx, center in enumerate(centers):
		list(center)
		if np.nan in center:
			centers[centerIdx] = str('NaN')
		else:
			center_ints = [int(round(dist)) for dist in center]
			centers[centerIdx] = center_ints		
else:
	centers = center_of_masss(mask)
	list(centers)	
	enters = [int(round(dist)) for dist in centers]
'''

out_img = nib.Nifti1Image(mask4D, img.affine, img.header)
out_dir = re.sub("native2std.*", "", input_img)
out_path = out_dir+'4dClusters.nii.gz'
nib.save(out_img, out_path)
	
print("Cluster_Index\tN_Voxels\tCenter_Coordinate")
for idx,cluster in enumerate(cleaned_clusters):
	print(str(idx)+"\t"+str(len(cluster))+"\t"+str(centers[idx]))

