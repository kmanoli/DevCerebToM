#!/usr/bin/env python3

import os
import glob
import pandas as pd
import nibabel as nib
import numpy as np
from nilearn.glm.second_level import SecondLevelModel
from nilearn.glm import threshold_stats_img

####################
### ALL CHILDREN
####################

# Set data directory 
data_dir = r'/data/project/DevCerebToM/children/'

# Create a directory for results
out_dir = os.path.join(data_dir, 'secondlevel_stv_results')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
    
# Define seed
seed = "r_crus_I_all"

# Filenames for each subject's seed-to-voxel maps for each seed
fnames = glob.glob(os.path.join(data_dir, 'firstlevel_stv_results', '**', f'{seed}*unthresh.nii.gz'), recursive=True)

# Set number of subjects and specify design matrix (one-sample t-test so one column for intercept only)
n_subs = len(fnames)
design_matrix = pd.DataFrame([1] * n_subs, columns=['intercept'])

# Fit second level model
second_level_model = SecondLevelModel().fit(
        fnames, design_matrix=design_matrix)
    
# Compute the contrast (one-sample t-test)
z_map = second_level_model.compute_contrast(output_type='z_score')

    # FDR thresholded map
thresholded_map, threshold = threshold_stats_img(
        z_map, alpha=.05, height_control='fdr')

# Save 
nib.save(z_map, os.path.join(out_dir, f"all_children_{seed}_stv_one_sample_unthresh.nii"))
nib.save(thresholded_map, os.path.join(out_dir, f"all_children_{seed}_stv_one_sample_fdr.nii"))

####################
### ToM PASSERS
####################

# Define subjects
sub_list_path = os.path.join(data_dir, 'passers.txt')  
with open(sub_list_path, 'r') as file:
    subjects = [line.strip() for line in file]

# Define seed
seed = "r_crus_II_pass"

for sub in subjects:  # Iterate over passers

        # Initialize an empty list to store filenames for selected subjects
        selected_fnames = []

        # Get data files (search every subject's subdirectory)
        fnames = glob.glob(os.path.join(data_dir, 'firstlevel_stv_results', sub, f'{seed}*unthresh.nii.gz'))

        # Extend selected_fnames with fnames
        selected_fnames.extend(fnames)

# Set number of subjects and specify design matrix (one-sample t-test so one column for intercept only)
n_subs = len(selected_fnames)
design_matrix = pd.DataFrame([1] * n_subs, columns=['intercept'])

# Fit second level model
second_level_model = SecondLevelModel().fit(
        selected_fnames, design_matrix=design_matrix)
    
# Compute the contrast (one-sample t-test)
z_map = second_level_model.compute_contrast(output_type='z_score')

    # FDR thresholded map
thresholded_map, threshold = threshold_stats_img(
        z_map, alpha=.05, height_control='fdr')

# Create a directory for results
out_dir_pass = os.path.join(data_dir, 'pass', 'secondlevel_stv_results')
if not os.path.exists(out_dir):
     os.makedirs(out_dir)

# Save 
nib.save(z_map, os.path.join(out_dir_pass, f"pass_{seed}_stv_one_sample_unthresh.nii"))
nib.save(thresholded_map, os.path.join(out_dir_pass, f"pass_{seed}_stv_one_sample_fdr.nii"))