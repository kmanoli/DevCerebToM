# Functional recruitment and connectivity of the cerebellum supports the emergence of Theory of Mind in early childhood

## This repository contains scripts and data to reproduce analyses preprinted here: https://www.biorxiv.org/content/10.1101/2024.04.02.586955v1

## Instructions:
If you'd like to reproduce the presented analyses, you can run the following scripts:

### A. Contrast analyses

- **1. /scripts/contrasts/<sample_name>/firstlevel_contrast.py:** These scripts generate first-level contrast maps for the contrast ToM > pain for each sample (Richardson adults, CCC adults, and children).
    _Expected outputs_: First-level contrast maps for each subject
      _Expected runtime_: ~1 day per sample

- **2. /scripts/cereb_isolation/cereb_isolation.m:** This script isolates the anatomical and first-level contrast functional images from the cerebral cortex and normalizes them to SUIT space for each sample.
    _Expected outputs_: Isolated and normalized cerebellum from T1-weighted and first-level contrast functional images for each subject
      _Expected runtime_: ~1 day per sample

- **3. /scripts/contrasts/<sample_name>/secondlevel_<analysis_type>.py:** These scripts generate second-level maps for each sample (Richardson adults, CCC adults, and children). We include sripts for group-level one-sample t-tests in each adult sample, the entire developmental sample, as well as children with and without ToM abilities. We provide an additional script for a group-level GLM in the entire developmental sample, in which children's ToM score is added as a covariate.
    _Expected outputs_: Second-level one-sample t-tests for each sample ('results_and_figures/contrast_flatmaps/<sample_name>/one_sample_fdr.nii') or GLM ('results_and_figures/contrast_flatmaps/children/GLM_fdr.nii')
      _Expected runtime_: ~5-10 minutes per sample

<img width="428" alt="Screen Shot 2024-05-14 at 4 25 14 pm" src="https://github.com/kmanoli/DevCerebToM/assets/44278225/d66817e5-483a-4652-8990-808d7023296c">

### B. Seed-to-voxel functional connectivity

- **1. /scripts/seed_to_voxel/first_level/<sample_name>_firstlevel_stv.py:** These scripts generate first-level seed-to-voxel functional connectivity maps for Richardson adults and children.
    _Expected outputs_: First-level seed-to-voxel functional connectivity maps for each subject
      _Expected runtime_: ~1 day per sample

- **2. /scripts/seed_to_voxel/one_sample/one_sample_stv.py:** This script generates second-level seed-to-voxel functional connectivity maps for children (one-sample t-tests).
    _Expected outputs_: Second-level one-sample t-test maps for children (Outputs for each seed in: 'results_and_figures/seed_to_voxel/<seed_name>')
      _Expected runtime_: ~15-30 minutes

- **3. /scripts/seed_to_voxel/two_sample/two_sample_stv.py:** This script generates second-level seed-to-voxel functional connectivity maps for children with vs. withput ToM abilities, and children vs. adults (two-sample t-tests).
    _Expected outputs_: Second-level two-sample t-test maps for adults and children (Outputs for each seed in: 'results_and_figures/seed_to_voxel/<seed_name>')
      _Expected runtime_: ~15-30 minutes

<img width="461" alt="Screen Shot 2024-05-14 at 4 35 56 pm" src="https://github.com/kmanoli/DevCerebToM/assets/44278225/add251a3-f8ce-411b-b889-3deb7f9da5ac">

### C. Dynamic causal modelling

- **1. /scripts/DCM/ROI_definition/<sample_name>/SPM_firstlevel.m:** These scripts generate first-level contrast maps for the contrast ToM > pain for each sample (Richardson adults, CCC adults, and children). This is a precursor for ROI definition and DCM analyses (as we need the SPM_results.mat file for each subject). In the case of CCC adults, we have included a script to smooth the data for consistency with the Richardson dataset ('/scripts/DCM/ROI_definition/adults_ccc/smooth.m')
    _Expected outputs_: First-level contrast maps for each subject
      _Expected runtime_: ~1 day per sample

- **2. /scripts/DCM/ROI_definition/<sample_name>/roi_def.m:** These scripts generate ROIs for DCM analyses. Local maxima and further instructions for ROI definition are provided within the scripts.
    _Expected outputs_: ToM ROIs for each subject
      _Expected runtime_: ~1 day per sample

- **3. /scripts/DCM/run_DCM/run_DCM.m:** This script runs first-level DCM analysis. Further instructions are provided within the script.
    _Expected outputs_: First-level DCM matrices for each subject
      _Expected runtime_: ~3-4 days per sample

- **4. /scripts/DCM/run_DCM/run_peb_2ndlevel.m:** This script runs second-level DCM analysis (Parametric Empirical Bayes model). Further instructions are provided within the script.
    _Expected outputs_: Second-level DCM matrices for each sample (Output in tabular format in: 'results_and_figures/DCM_heatmaps/<sample_name>')
      _Expected runtime_: ~3-4 days per sample

<img width="652" alt="Screen Shot 2024-05-14 at 4 47 44 pm" src="https://github.com/kmanoli/DevCerebToM/assets/44278225/29d62caa-24d7-4ad1-a1d0-a24b6e63d877">

### D. Figures

- **1. /results_and_figures/contrast_flatmaps.ipynb:** This Jupyter notebook generates cerebellar flatmaps for second-level contrast analyses. 
      _Expected runtime_: ~10 minutes

- **2. /results_and_figures/seed_to_voxel_figs.ipynb:** This Jupyter notebook generates cerebral cortex maps for second-level seed-to-voxel functional connectivity. 
      _Expected runtime_: ~15 minutes

- **3. /results_and_figures/DCM_heatmaps.ipynb:** This Jupyter notebook generates heatmaps for second-level DCM results.
      _Expected runtime_: ~10 minutes 


## Data:
This study is based on openly accessible data hosted on OpenNeuro. 

You can download the data here: 

- Richardson et al., 2018: https://openneuro.org/datasets/ds000228/versions/1.1.0
- CCC: https://openneuro.org/datasets/ds003798


## Requirements:
Matlab scripts were run in Matlab 2022b, Python scripts were run in Python 3.9.7 on Linux.

_Matlab_: 
Our analysis code makes use of open software. DCM analysis (including contrast analyses and ROI definition) was run in SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/). Cerebellar isolation and normalization was run in SUIT (v. 3.5; https://github.com/jdiedrichsen/suit/releases/tag/3.5)(https://github.com/kmanoli/DevCerebToM/assets/44278225/0afa6301-2de1-4923-9a49-1c1bca37a8be). 

_Python_:
We made use of the following packages: nilearn 0.10.1, nibabel 5.1.0, pandas 2.0.3, numpy 1.25.2, SUITpy 1.2.0, seaborn 0.12.2, and matplotlib 3.7.2. Surface visualization was done with surfplot 0.2.0 and neuromaps 0.0.4, after projecting volumetric images to cortical surfaces using Workbench Command (https://www.humanconnectome.org/software/connectome-workbench). 

The scripts can be run on a standard desktop computer, however, we recommend running these analyses (especially first-level models and DCM analyses) on a cluster to parallelize computations.

## Install: 
Follow these instructions to install Matlab (license required): https://de.mathworks.com/help/install/, Python: https://www.python.org/downloads/, and Workbench Command (for surface visualizations): https://www.humanconnectome.org/software/get-connectome-workbench. Then, clone this repository:
git clone https://github.com/kmanoli/DevCerebToM.git
(This should take a few seconds only.)


