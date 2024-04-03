
% This script assumes that you have specified SPM results .mat files for
% the contrast ToM > pain (adjusted for effects of interest) - see
% SPM_firstlevel.m script. 

% Here, we define ROIs with a progressively tolerant threshold. Subject
% ROI definition was manually inspected. As a tip, if SPM gives a warning
% during ROI creation for a subject (output in orange), then a different
% threshold should be used for this subject (e.g., .10 instead of .05, or 1
% instead of .10).


% ROI NAMES AND COORDINATES FOR INPUT BELOW

% pc - [0 – 60 40]
% dmpfc - [0 50 35]
% vmpfc - [0 50 5]
% rtpj - [50 – 55 25]
% ltpj - [-50 – 55 25]
% rcus2 - [26 -80 -39]
% lcrus2 - [-26 -78 -38]

%--------------------------------------------------------------------------
% Initialize SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% Specify data path
data_dir = '/data/project/DevCerebToM/DCM_dirs/adults_ccc';

% Subject list
subject_list_file = fullfile(data_dir, 'adults_ccc.txt');
subjects = importdata(subject_list_file);

% Specify thresholded SPM and extract VOI data for each subject

for subject = subjects
    
subject = num2str(subject, '%03d'); % Zero-pads each number so that the subject ID is 3 characters long


% Insert the subject's SPM .mat filename here
sub_dir = fullfile(data_dir, ['data/sub_' subject]);
spm_mat_file = fullfile(sub_dir,  'SPM_results/SPM.mat')

% Start batch
clear matlabbatch;
matlabbatch{1}.spm.util.voi.spmmat  = cellstr(spm_mat_file);
matlabbatch{1}.spm.util.voi.adjust  = 2;                    % Effects of interest contrast number (position in SPM.mat)
matlabbatch{1}.spm.util.voi.session = 1;                    % Session index
matlabbatch{1}.spm.util.voi.name    = 'lcrus2';             % VOI name

% Define thresholded SPM for finding the subject's local peak response
matlabbatch{1}.spm.util.voi.roi{1}.spm.spmmat      = {''};
matlabbatch{1}.spm.util.voi.roi{1}.spm.contrast    = 1;     % Index of contrast for choosing voxels (position of ToM-pain contrast in SMP.mat)
matlabbatch{1}.spm.util.voi.roi{1}.spm.conjunction = 1;
matlabbatch{1}.spm.util.voi.roi{1}.spm.threshdesc  = 'none';
matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh      = 0.05; % alpha level

% For subjects where alpha=.05 returns no voxels/has an incorrect center
%matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh      = 0.1; % alpha level
% For subjects where alpha=.10 returns no voxels/has an incorrect center
%matlabbatch{1}.spm.util.voi.roi{1}.spm.thresh      = 1; % alpha level

matlabbatch{1}.spm.util.voi.roi{1}.spm.extent      = 5;
matlabbatch{1}.spm.util.voi.roi{1}.spm.mask ...
    = struct('contrast', {}, 'thresh', {}, 'mtype', {});

% Define large fixed outer sphere (meta-analytic group VOI)
matlabbatch{1}.spm.util.voi.roi{2}.sphere.centre     = [-24 -78 -38]; % Set coordinates here

% CORTEX
%matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius     = 15;            % Radius (mm)
% CEREBELLUM
matlabbatch{1}.spm.util.voi.roi{2}.sphere.radius     = 11;            % Radius (mm)

matlabbatch{1}.spm.util.voi.roi{2}.sphere.move.fixed = 1;

% Define smaller inner sphere which jumps to the peak of the outer sphere
% (subject-specific VOI)
matlabbatch{1}.spm.util.voi.roi{3}.sphere.centre           = [0 0 0]; % Leave this at zero since centre will move freely

% CORTEX
%matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius           = 8;       % Set radius here (mm)
% CEREBELLUM
matlabbatch{1}.spm.util.voi.roi{3}.sphere.radius           = 5;       % Set radius here (mm)

matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.spm  = 1;       % Index of SPM within the batch (here:1 because thresholded SPM was defined before spheres)
matlabbatch{1}.spm.util.voi.roi{3}.sphere.move.global.mask = 'i2';    % Index of the outer sphere within the batch

% Include voxels in the thresholded SPM (i1) and the mobile inner sphere (i3)
matlabbatch{1}.spm.util.voi.expression = 'i1 & i3'; 

% Run the batch
spm_jobman('run',matlabbatch);
end 