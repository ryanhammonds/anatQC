#!/usr/bin/env python

import sys
import re
import numpy as np
import nibabel as nib
from scipy.ndimage import center_of_mass

# Gets nifti as input
input_img = sys.argv[1]
voxel_std = float(sys.argv[2])
clust_size = int(sys.argv[3])

# Variables for naming saved outputs
out_dir = re.sub("native2std.*", "", input_img)
subj_id = re.sub(".*sub", "sub", input_img)
subj_id = re.sub("/.*", "", subj_id)
base_dir = re.sub("sub.*", "", input_img)

# Load nifti to object
img = nib.load(input_img)
img_data = img.get_data()
data_obj = img.dataobj

# Get dimensions
x_len = data_obj.shape[0]
y_len = data_obj.shape[1]
z_len = data_obj.shape[2]

# Creates a list of all nonzero intensities, used to calculate std/mean
voxelValsList = []
for z in range(0, z_len):
    for y in range(0, y_len):
        for x in range(0, x_len):
            if img_data[x, y, z] > 0:
                voxelValsList.append(img_data[x, y, z])

meanIntensity = np.mean(voxelValsList)
stdIntensity = np.std(voxelValsList)
threshold = meanIntensity - (voxel_std * stdIntensity)

# For std arguments that are too exclusive (negative threshold value)
while threshold < 0 and voxel_std >= 1:
    voxel_std = voxel_std - 0.1
    threshold = meanIntensity - (voxel_std * stdIntensity)

voxel_write = out_dir+'voxel_thr_used.txt'
print(voxel_std, file=open(voxel_write, "w"))

# List of all coordinates below threshold
flagList = []
for z in range(0, z_len):
    for y in range(0, y_len):
        for x in range(0, x_len):
            if img_data[x, y, z]:
                if img_data[x, y, z] < threshold:
                    flagList.append([x, y, z])

# Recursive function, finds clusters in 3D
sys.setrecursionlimit(len(flagList))


def checkContig(pt, flagList, clustList):
    x = pt[0]
    y = pt[1]
    z = pt[2]
    if pt in flagList:
        clustList.append(pt)
        flagList.remove(pt)

    # All voxels surrounding a given voxels (ie rubic's cube)
    points = [
        [x, y, z-1], [x+1, y, z-1], [x+1, y+1, z-1], [x+1, y-1, z-1],
        [x-1, y, z-1], [x-1, y+1, z-1], [x-1, y-1, z-1], [x, y+1, z-1],
        [x, y-1, z-1], [x, y, z+1], [x+1, y, z+1], [x+1, y+1, z+1],
        [x+1, y-1, z+1], [x-1, y, z+1], [x-1, y+1, z+1], [x-1, y-1, z+1],
        [x, y+1, z+1], [x, y-1, z+1], [x, y-1, z], [x+1, y-1, z],
        [x+1, y-1, z+1], [x+1, y-1, z-1], [x-1, y-1, z], [x-1, y-1, z+1],
        [x-1, y-1, z-1], [x, y-1, z+1], [x, y-1, z-1], [x, y+1, z],
        [x+1, y+1, z], [x+1, y+1, z+1], [x+1, y+1, z-1], [x-1, y+1, z],
        [x-1, y+1, z+1], [x-1, y+1, z-1], [x, y+1, z+1], [x, y+1, z-1],
        [x-1, y, z], [x-1, y+1, z], [x-1, x+1, z+1], [x-1, y+1, z-1],
        [x-1, y-1, z], [x-1, y-1, z+1], [x-1, y-1, z-1], [x-1, y, z+1],
        [x-1, y, z-1], [x+1, y, z], [x+1, y+1, z], [x+1, y+1, z+1],
        [x+1, y+1, z-1], [x+1, y-1, z], [x+1, y-1, z+1], [x+1, y-1, z-1],
        [x+1, y, z+1], [x+1, y, z-1]
    ]

    # Remove coords adjacent to a voxel w/intensity of 0.
    for coord in points:
        if img_data[coord[0], coord[1], coord[2]] == 0:
            clustList.remove(pt)
            break

    # Recusion with adjacent points
    for coord in points:
        if coord in flagList and coord not in clustList:
            pt = list(coord)
            clustList.append(pt)
            flagList.remove(pt)
            checkContig(pt, flagList, clustList)
    return clustList


# Flag clusters
clusters = []
for xyz in flagList:
    singleCluster = []
    singleCluster = checkContig(xyz, flagList, singleCluster)
    clusters.append(singleCluster)

# Cluster thresholding
cleaned_clusters = []
for cluster in clusters:
    if len(cluster) > clust_size:
        cleaned_clusters.append(cluster)

# Make mask and info text file if any cluster survives
if len(cleaned_clusters) > 0:
    # Create 4D mask of clusters and output sizes and centers
    centers = []
    n_voxels = []
    firstImg = True
    for idx, cluster in enumerate(cleaned_clusters):
        # Create mask
        mask = np.zeros(np.shape(img_data))
        for voxel in cluster:
            mask[voxel[0], voxel[1], voxel[2]] = 1

        # Find center of cluster, if center isn't in mask, delete cluster
        center = center_of_mass(mask)
        list(center)
        center = [int(round(dist)) for dist in center]
        if mask[center[0], center[1], center[2]] != 0:
            voxels = len(cleaned_clusters[idx])
            n_voxels.append(voxels)
            centers.append(center)
            mask = mask[..., np.newaxis]
            if firstImg is True:
                mask4D = mask[..., np.newaxis]
            elif firstImg is False:
                mask4D = np.concatenate((mask4D, mask), axis=3)
            firstImg = False

    # Save nifti
    if len(centers) > 0:
        out_img = nib.Nifti1Image(mask4D, img.affine, img.header)
        out_path = out_dir+'4dClusters.nii.gz'
        nib.save(out_img, out_path)
        # Save .tsv file
        voxel_str = str(voxel_std)
        clust_str = str(clust_size)
        tsv = out_dir+'4dClusters_'+'vox'+voxel_str+'_clust'+clust_str+'.tsv'
        header = str("Cluster_Index\tN_Voxels\tCenter_Coordinate")
        print(header, file=open(tsv, "w"))
        for idx, center in enumerate(centers):
            print(str(idx)+"\t"+str(n_voxels[idx])+"\t"+str(center), file=open(tsv, "a"))
