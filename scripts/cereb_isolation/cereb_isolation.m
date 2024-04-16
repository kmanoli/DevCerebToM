% Before running this script, make sure to initialise SPM by running 'spm
% fmri' (without quotes) in the command line

% Work in a clean MATLAB workspace
clearvars;

% Add input folder to path
addpath(genpath('/data/project/DevCerebToM/suit/input/'));

% Make subject list
subjects = importdata('/data/project/DevCerebToM/suit/input/sub_list.txt')
subjects = cellstr(num2str(subjects));
subjects = subjects';

%% ISOLATION/SEGMENTATION

for subject = 1:length(subjects)

% Check whether the (unsmoothed) T1w images have been unzipped. If not, unzip them using
% gunzip
if isfile(['/data/project/DevCerebToM/suit/input/' subjects{subject} '/' subjects{subject} '_T1w.nii']) == 0
    display('Image has not been unzipped; unzipping now')
    gunzip(['/data/project/DevCerebToM/suit/input/' subjects{subject} '/' subjects{subject} '_T1w.nii.gz'])
else
    display('Image is already unzipped; doing nothing')
end


% Isolate and segment the cerebellum
suit_isolate_seg({['/data/project/DevCerebToM/suit/input/' subjects{subject} '/' subjects{subject} '_T1w.nii']})

% Isolation mask of the cerebellum (output file ending in '_pcereb.nii')
% should be inspected visually and hand-corrected (e.g., in MRICroN)
end 

%% NORMALIZATION

% Normalize the isolated cerebellum to SUIT atlas
for subject = 1:length(subjects)
    job.subjND(subject).gray = {['/data/project/DevCerebToM/suit/input/' subjects{subject} '/' subjects{subject} '_T1w_seg1.nii']}; 
    job.subjND(subject).white = {['/data/project/DevCerebToM/suit/input/' subjects{subject} '/' subjects{subject} '_T1w_seg2.nii']};
    job.subjND(subject).isolation = {['/data/project/DevCerebToM/suit/input/' subjects{subject} '/c_' subjects{subject} '_T1w_pcereb.nii']};
end

suit_normalize_dartel(job);

%% RESLICING

% Reslice cerebellum to SUIT space
for subject=1:length(subjects)
    job.subj(subject).affineTr = {['/data/project/DevCerebToM/suit/input/' subjects{subject} '/Affine_' subjects{subject} '_T1w_seg1.mat']}; 
    job.subj(subject).flowfield = {['/data/project/DevCerebToM/suit/input/' subjects{subject} '/u_a_' subjects{subject} '_T1w_seg1.nii']};
    job.subj(subject).resample = {['/data/project/DevCerebToM/suit/input/' subjects{subject} '/' subjects{subject} '_T1w_.nii']};
    job.subj(subject).mask = {['/data/project/DevCerebToM/suit/input/' subjects{subject} '/c_' subjects{subject} '_T1w_pcereb.nii']};
end

suit_reslice_dartel(job);

