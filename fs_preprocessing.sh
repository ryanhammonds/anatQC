#!/bin/bash

# Extracts volumes and thicknessess from freesurfer recons and formats to csv.

set -e
# Arguments
recon_base=$1
output=$2

# Setup
if [[ ! -d $recon_base ]]; then
  echo 'Recon directory does not exist.'
  exit 1
fi

# Remove trailing slash.
recon_base=${recon_base%/}
output=${output%/}

if [[ ! -d $output ]]; then
  mkdir $output
fi

# Extract stats function (within subjects)
extractStats(){
  [[ ! -d $output/$1 ]] && mkdir $output/$1
  cat $recon_base/$1/stats/aseg.stats | sed "/#/d" | awk '{print $4}' > $output/$1/aseg_vol.txt
  for hemi in lh rh; do
    cat $recon_base/$1/stats/$hemi.aparc.stats | sed "/#/d" | awk '{print $4}' > $output/$1/$hemi.aparc_vol.txt
    cat $recon_base/$1/stats/$hemi.aparc.stats | sed "/#/d" | awk '{print $5}' > $output/$1/$hemi.aparc_thick.txt
  done
  cat $recon_base/$1/stats/aseg.stats | grep 'Estimated Total Intracranial Volume' | sed "s/, mm^3.*//" | sed "s/.*. //" > $output/$1/ICV.txt

  cat $output/$1/aseg_vol.txt $output/$1/lh.aparc_vol.txt $output/$1/rh.aparc_vol.txt > $output/$1/volumes.txt
  cat $output/$1/lh.aparc_thick.txt $output/$1/rh.aparc_thick.txt > $output/$1/thickness.txt
  cat $output/$1/volumes.txt $output/$1/thickness.txt > $output/$1/volumes+thickness.txt
}

# Create array for subjects who have all stats files
subjArray=()
for recon in $recon_base/*; do
  stats=$recon/stats
  if [[ -f $stats/aseg.stats && -f $stats/lh.aparc.stats && -f $stats/rh.aparc.stats ]]; then
    subj=$(echo $recon | sed "s#.*/##")
    subjArray+=($subj)
  fi
done

# Call extractStats function serially
for subj in ${subjArray[@]}; do
  extractStats $subj
done

echo -e "Rows = Features\nColumns = Subjects" >> $output/README.txt

# Create headers
cat $recon_base/${subjArray[0]}/stats/aseg.stats | sed "/#/d" | awk '{print $5}' > $output/aseg_labels.txt
cat $recon_base/${subjArray[0]}/stats/lh.aparc.stats | sed "/#/d" | awk '{print $1}' | sed -e "s/^/lh_/" > $output/lh.aparc_labels.txt
cat $recon_base/${subjArray[0]}/stats/rh.aparc.stats | sed "/#/d" | awk '{print $1}' | sed -e "s/^/rh_/"> $output/rh.aparc_labels.txt
cat $output/aseg_labels.txt $output/lh.aparc_labels.txt $output/rh.aparc_labels.txt > $output/volumes_header.txt
cat $output/lh.aparc_labels.txt $output/rh.aparc_labels.txt > $output/thickness_header.txt
cat $output/volumes_header.txt $output/thickness_header.txt > $output/volumes+thickness_header.txt
rm $output/aseg_labels.txt $output/lh.aparc_labels.txt $output/rh.aparc_labels.txt

# Merge volumes and thicknesses (between subjects)
paste -d',' $output/*/volumes.txt > $output/volumes.csv
paste -d',' $output/*/thickness.txt > $output/thickness.csv
paste -d',' $output/*/volumes+thickness.txt > $output/volumes+thickness.csv
paste -d',' $output/*/ICV.txt >> $output/ICV.txt
