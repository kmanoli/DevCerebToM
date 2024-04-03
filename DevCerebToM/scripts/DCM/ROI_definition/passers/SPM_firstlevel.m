% Initialize SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% Specify data path
data_dir = '/data/project/DevCerebToM/DCM_dirs/nonpassers';

% Subject list
subject_list_file = fullfile(data_dir, 'nonpassers.txt');
subjects = importdata(subject_list_file);

% Specify and estimate 1st level model and contrasts

for subject = subjects
    
subject = num2str(subject, '%03d'); % Zero-pads each number so that the subject ID is 3 characters long

% Unzip files
files = dir(fullfile(data_dir, ['data/sub_' subject], '*bold.nii.gz'));
if ~isempty(files)
    % Unzip files if they exist
    for i = 1:numel(files)
        zip_file = fullfile(data_dir, ['data/sub_' subject], files(i).name); % Corrected variable name
        fprintf('Image has not been unzipped; unzipping %s\n', zip_file);
        gunzip(zip_file);
    end
else
    fprintf('No images to unzip for subject %s\n', subject);
end

% Turn nuisance regressors into .txt files
conf = dir(fullfile(data_dir, ['data/sub_' subject], '*ART_and_CompCor_nuisance_regressors.mat'));

if ~isempty(conf)
    % Assuming there's only one matching file, you can directly access the first file
    file_path = fullfile(data_dir, ['data/sub_' subject], conf(1).name);
    load(file_path);
else
    fprintf('No matching file found for subject %s\n', subject);
end

writematrix(R, fullfile(data_dir, ['data/sub_' subject], 'nuisance_regressors.txt'));

%%%%%%%%%%

% Specify data path
sub_dir = fullfile(data_dir, ['data/sub_' subject]);

% Specify scans and volumes
scans = spm_select('ExtFPList', fullfile(data_dir) ,'.*.nii$', 1:168);

clear matlabbatch

% Model preliminaries
matlabbatch{1}.spm.stats.fmri_spec.dir = {fullfile(sub_dir,  'SPM_results')};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

% Select scans for model specification
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(scans); % Pass scans as cell

% Specify condition onsets and durations
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name = 'tom';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset = [86
                                                         98
                                                         120
                                                         176
                                                         238
                                                         252
                                                         300];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration = [4
                                                            6
                                                            4
                                                            16
                                                            6
                                                            8
                                                            6];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name = 'pain';
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset = [70
                                                         92
                                                         106
                                                         136
                                                         194
                                                         210
                                                         228
                                                         262
                                                         312];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration = [4
                                                            2
                                                            4
                                                            10
                                                            4
                                                            12
                                                            6
                                                            6
                                                            4];
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod = 0;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod = struct('name', {}, 'param', {}, 'poly', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});

% Add nuisance regressors
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {fullfile(sub_dir, 'nuisance_regressors.txt')};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

% Model estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Contrast estimation
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'tom-pain';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.fcon.name = 'effects of interest';
matlabbatch{3}.spm.stats.con.consess{2}.fcon.weights = [1 0
                                                        0 1];
matlabbatch{3}.spm.stats.con.consess{2}.fcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

% Run job
spm_jobman('run', matlabbatch);

end