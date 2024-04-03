#!/usr/bin/env python3

import os
import glob
import nibabel as nib
import pandas as pd 
import numpy as np 
from nilearn.maskers import NiftiMasker
from nilearn.maskers import NiftiSpheresMasker
from pymatreader import read_mat

# Set directory for data loading
data_dir = r'/data/project/DevCerebToM/adults_richardson/'

# Define subjects
sub_list = os.path.join(data_dir, 'adults_richardson.txt')  

with open(sub_list, 'r') as file:
    subjects = [line.strip() for line in file]

# Create a directory for results
out_dir = os.path.join(data_dir, 'firstlevel_stv_results')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Run first-level model 
for sub in subjects:

    # Make results subdirectory for every subject
    out_subject_dir = os.path.join(out_dir, sub)
    if not os.path.exists(out_subject_dir):
        os.mkdir(out_subject_dir)

     # Load image for single subject
    scan_files = glob.glob(os.path.join(data_dir, 'data', sub, '*.nii.gz'))
    if scan_files:
        scan_file = scan_files[0]
        print(f"Running first-level model for subject {sub}")
        scan = nib.load(scan_file)
    
        # Load nuisance regressors
        confound_file = glob.glob(os.path.join(data_dir, f'data/{sub}/conf/*ART_and_CompCor_nuisance_regressors.mat'))[0]
        conf_mat = read_mat(confound_file)
        confounds = pd.DataFrame(conf_mat['R'])

        # Define seeds
        seed_list = [
            {"name": "r_crus_II",
            "seed": [(26,-80,-39)]},
            {"name": "l_crus_II",
            "seed": [(-26,-78,-38)]},
        ]

        for seed in seed_list:
            # Extract timeseries for seeds 
            seed_masker = NiftiSpheresMasker(
                seed["seed"],
                radius=5,
                detrend=True,
                standardize="zscore_sample",
                standardize_confounds="zscore_sample",
                low_pass=0.1,
                high_pass=0.01,
                t_r=2,
                memory="nilearn_cache",
                memory_level=1,
                verbose=0,
            )
            seed_time_series = seed_masker.fit_transform(
            scan, confounds=[confounds]
            )

            # Extract timeseries for rest of the brain 
            brain_masker = NiftiMasker(
            detrend=True,
            standardize="zscore_sample",
            standardize_confounds="zscore_sample",
            low_pass=0.1,
            high_pass=0.01,
            t_r=2,
            memory="nilearn_cache",
            memory_level=1,
            verbose=0,
            )

            brain_time_series = brain_masker.fit_transform(
            scan, confounds=[confounds]
            )

            # Perform seed-to-voxel correlation
            seed_to_voxel_correlations = (
            np.dot(brain_time_series.T, seed_time_series) / seed_time_series.shape[0]
            )

            # Turn correltation map to NIfTI and save
            seed_to_voxel_correlations_img = brain_masker.inverse_transform(
            seed_to_voxel_correlations.T
            )

            nib.save(seed_to_voxel_correlations_img, os.path.join(out_dir, f"{sub}/{seed['name']}_stv_unthresh.nii.gz"))

        # Track execution
        print(f"Finished running first-level model for subject {sub}")
    else:
        print(f"Subject {sub}: Scan not found") 