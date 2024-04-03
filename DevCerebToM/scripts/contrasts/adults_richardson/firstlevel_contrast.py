#!/usr/bin/env python3

# Imports
import os
import glob
import pandas as pd
import nibabel as nib
from numpy import array
from nilearn.glm.first_level import FirstLevelModel
from pymatreader import read_mat

# Set data directory 
data_dir = r'/data/project/DevCerebToM/adults_richardson'

# Define subjects
sub_list = os.path.join(data_dir, 'adults_richardson.txt')  

with open(sub_list, 'r') as file:
    subjects = [line.strip() for line in file]

# Create a directory for results
out_dir = os.path.join(data_dir, 'firstlevel_results')
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# Run first-level model 
for sub in subjects:

    # Load image for single subject
    scan_files = glob.glob(os.path.join(data_dir, 'data', sub, '*.nii.gz'))
    if scan_files:
        scan_file = scan_files[0]
        print(f"Running first-level model for subject {sub}")
        scan = nib.load(scan_file)

        # Load events
        events = pd.read_table(os.path.join(data_dir, 'movie_events_seconds.tsv'))

        # Load nuisance regressors
        confound_file = glob.glob(os.path.join(data_dir, f'data/{sub}/conf/*ART_and_CompCor_nuisance_regressors.mat'))[0]
        conf_mat = read_mat(confound_file)
        confounds = pd.DataFrame(conf_mat['R'])

        # Define first level model, apply it to scan, and specify the design matrix
        fmri_glm = FirstLevelModel(t_r=2, standardize=False, hrf_model='spm + derivative')
        fmri_glm = fmri_glm.fit(scan, events, confounds=confounds)
        design_matrix = fmri_glm.design_matrices_[0]

        # Specify conditions and calculate the ToM-pain contrast
        pain = array([0] * len(design_matrix.axes[1]))
        pain[0] = 1
        tom = array([0] * len(design_matrix.axes[1]))
        tom[2] = 1

        conditions = {
            'pain': pain,
            'tom':  tom,
        }

        tom_minus_pain = conditions['tom'] - conditions['pain']

        z_map = fmri_glm.compute_contrast(tom_minus_pain, output_type='z_score')

        # Make results subdirectory for every subject
        out_subject_dir = os.path.join(out_dir, sub)
        if not os.path.exists(out_subject_dir):
            os.mkdir(out_subject_dir)

        # Save unthresholded map
        nib.save(z_map, os.path.join(out_dir, f"{sub}/{sub}_tom_vs_pain_unthresh.nii.gz"))

        # Track execution
        print(f"Finished running first-level model for subject {sub}")
    else:
        print(f"Subject {sub}: Scan not found")

    
