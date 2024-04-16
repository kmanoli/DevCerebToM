#!/usr/bin/env python3

# Imports
import os
import glob
import pandas as pd
import nibabel as nib
from nilearn.glm.second_level import SecondLevelModel
from nilearn.glm import threshold_stats_img

# Set data directory 
data_dir = r'/data/project/DevCerebToM/adults_richardson/'

# List directories for cortex and cerebellum
# Note: This assumes that you have isolated the first-level images with the script 'cereb_isolation.m'
# and have placed them in a 'firstlevel_results_cer' subfolder
dir_list = [os.path.join(data_dir, 'firstlevel_results'), 
            os.path.join(data_dir, 'firstlevel_results_cer')]

# Run second-level model separately for the cerebellum and the cortex
for dir in dir_list:

    # Get data files (search every subject's subdirectory)
    fnames = glob.glob(os.path.join(data_dir, dir, '**', '*unthresh.nii.gz'), recursive=True)

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
    
    # Create a directory for results
    if dir == os.path.join(data_dir, 'firstlevel_results'):
        out_dir = os.path.join(data_dir, 'secondlevel_results')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    elif dir == os.path.join(data_dir, 'firstlevel_results_cer'):
        out_dir = os.path.join(data_dir, 'secondlevel_results_cer')
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    # Save unthresholded and thresholded maps
    nib.save(z_map, os.path.join(out_dir, 'one_sample_unthresh.nii.gz'))
    nib.save(thresholded_map, os.path.join(out_dir, 'one_sample_fdr.nii.gz'))