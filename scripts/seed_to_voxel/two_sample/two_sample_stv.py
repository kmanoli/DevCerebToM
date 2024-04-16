#!/usr/bin/env python3

import os
import glob
import pandas as pd
import nibabel as nib
import numpy as np
from nilearn.glm.second_level import SecondLevelModel
from nilearn.glm import threshold_stats_img

##################################
### ToM PASSERS VS. NON-PASSERS
##################################

# Set data directory 
data_dir = r'/data/project/DevCerebToM/children/'

# Create a directory for results
out_dir = os.path.join(data_dir, 'secondlevel_stv_results')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Import covariates
cov = pd.read_csv(os.path.join(data_dir, 'participants.tsv', sep='\t', usecols = ['FB_Composite', 'Child_Adult']))
cov = cov[cov['Child_Adult'] == 'child']

# Get values for FB task performance
tom_score = cov['FB_Composite'].values
tom_score = np.where(np.isin(tom_score, [0, 1, 2, 3]), 'fail', 'pass')

# Set number of subjects and specify design matrix 
n_subs = len(tom_score)
intercept = np.ones(n_subs)
design_matrix = pd.DataFrame(
    np.vstack((tom_score, intercept)).T,
    columns=['tom_score', 'intercept'],
)
    
# Create seed list
seed_list = ["r_crus_I_all", "r_crus_II_pass"]

for seed in seed_list: 
    # Filenames for each subject's seed-to-voxel maps for each ROI
    fnames = glob.glob(os.path.join(data_dir, 'firstlevel_stv_results', '**', f'{seed}*unthresh.nii'), recursive=True)
    fnames = sorted(fnames) # So that the covariate values are assigned to the correct subject

    # Fit second level model
    second_level_model = SecondLevelModel().fit(
            fnames, design_matrix=design_matrix)
        
    # Compute the contrast 
    z_map = second_level_model.compute_contrast(second_level_contrast='tom_score', output_type='z_score')
    
     # FDR thresholded map
    thresholded_map, threshold = threshold_stats_img(
         z_map, alpha=.05, height_control='fdr')

    # Save 
    nib.save(z_map, os.path.join(out_dir, f"{seed}_stv_two_sample_unthresh.nii.gz"))
    nib.save(thresholded_map, os.path.join(out_dir, f"{seed}_stv_two_sample_fdr.nii.gz"))


##################################
### ToM PASSERS VS. ADULTS
##################################
    
# Set data directory 
data_dir = r'/data/project/DevCerebToM/'

# Create covariates and binarize them
group_type = np.repeat(['pass', 'adults'], [22, 22])
group_type = np.where(group_type == "pass", -1, 1)

# Specify design matrix
n_subs = len(group_type)
intercept = np.ones(n_subs)
design_matrix = pd.DataFrame(
    np.vstack((group_type, intercept)).T,
    columns=['group_type', 'intercept'],
)
    
# Define seed
seed_list = "r_crus_II_adults"

# Filenames for each subject's seed-to-voxel maps for each ROI
# Define subjects
# Passers
sub_list_path = os.path.join(data_dir, 'children', 'passers.txt')  
with open(sub_list_path, 'r') as file:
    subjects = [line.strip() for line in file]

for sub in subjects:  # Iterate over passers

        # Initialize an empty list to store filenames for selected subjects
        fnames_p = []

        # Get data files (search every subject's subdirectory)
        fnames = glob.glob(os.path.join(data_dir, 'children', 'firstlevel_stv_results', sub, f'{seed}*unthresh.nii.gz'))

        # Extend selected_fnames with fnames
        fnames_p.extend(fnames)

# Adults
sub_list_path = os.path.join(data_dir, 'adults_richardson', 'adults_richardson.txt')  
with open(sub_list_path, 'r') as file:
    subjects = [line.strip() for line in file]

for sub in subjects:  # Iterate over passers

        # Initialize an empty list to store filenames for selected subjects
        fnames_a = []

        # Get data files (search every subject's subdirectory)
        fnames = glob.glob(os.path.join(data_dir, 'adults_richardson', 'firstlevel_stv_results', sub, f'{seed}*unthresh.nii.gz'))

        # Extend selected_fnames with fnames
        fnames_a.extend(fnames)

# Concatenate fnames
fnames = fnames_p + fnames_a

# Fit second level model
second_level_model = SecondLevelModel().fit(
        fnames, design_matrix=design_matrix)
    
# Compute the contrast (one-sample t-test)
z_map = second_level_model.compute_contrast(second_level_contrast='group_type', output_type='z_score')

    # FDR thresholded map
thresholded_map, threshold = threshold_stats_img(
        z_map, alpha=.05, height_control='fdr')

# Save 
nib.save(z_map, os.path.join(out_dir, f"{seed}_pass_stv_two_sample_unthresh.nii.gz"))
nib.save(thresholded_map, os.path.join(out_dir, f"{seed}_pass_stv_two_sample_fdr.nii.gz"))

##################################
### ToM NONPASSERS VS. ADULTS
##################################
    
# Set data directory 
data_dir = r'/data/project/DevCerebToM/'

# Create covariates and binarize them
group_type = np.repeat(['fail', 'adults'], [22, 22])
group_type = np.where(group_type == "fail", -1, 1)

# Specify design matrix
n_subs = len(group_type)
intercept = np.ones(n_subs)
design_matrix = pd.DataFrame(
    np.vstack((group_type, intercept)).T,
    columns=['group_type', 'intercept'],
)
    
# Define seed
seed_list = "r_crus_II_adults"

# Filenames for each subject's seed-to-voxel maps for each ROI
# Define subjects
# Passers
sub_list_path = os.path.join(data_dir, 'children', 'nonpassers.txt')  
with open(sub_list_path, 'r') as file:
    subjects = [line.strip() for line in file]

for sub in subjects:  # Iterate over nonpassers

        # Initialize an empty list to store filenames for selected subjects
        fnames_f = []

        # Get data files (search every subject's subdirectory)
        fnames = glob.glob(os.path.join(data_dir, 'children', 'firstlevel_stv_results', sub, f'{seed}*unthresh.nii.gz'))

        # Extend selected_fnames with fnames
        fnames_f.extend(fnames)

# Adults
sub_list_path = os.path.join(data_dir, 'adults_richardson', 'adults_richardson.txt')  
with open(sub_list_path, 'r') as file:
    subjects = [line.strip() for line in file]

for sub in subjects:  # Iterate over passers

        # Initialize an empty list to store filenames for selected subjects
        fnames_a = []

        # Get data files (search every subject's subdirectory)
        fnames = glob.glob(os.path.join(data_dir, 'adults_richardson', 'firstlevel_stv_results', sub, f'{seed}*unthresh.nii.gz'))

        # Extend selected_fnames with fnames
        fnames_a.extend(fnames)

# Concatenate fnames
fnames = fnames_f + fnames_a

# Fit second level model
second_level_model = SecondLevelModel().fit(
        fnames, design_matrix=design_matrix)
    
# Compute the contrast (one-sample t-test)
z_map = second_level_model.compute_contrast(second_level_contrast='group_type', output_type='z_score')

    # FDR thresholded map
thresholded_map, threshold = threshold_stats_img(
        z_map, alpha=.05, height_control='fdr')

# Save 
nib.save(z_map, os.path.join(out_dir, f"{seed}_fail_stv_two_sample_unthresh.nii.gz"))
nib.save(thresholded_map, os.path.join(out_dir, f"{seed}_fail_stv_two_sample_fdr.nii.gz"))



