% Initialize SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');


% Specify data path
data_dir = '/data/project/DevCerebToM/adults_ccc';

% Subject list
subject_list_file = fullfile(data_dir, 'adults_ccc.txt');
subjects = importdata(subject_list_file);

% Smooth images
for subject = subjects
    
subject = num2str(subject, '%03d'); % Zero-pads each number so that the subject ID is 3 characters long

% Unzip files
files = dir(fullfile(data_dir, ['sub_' subject], '*bold.nii.gz'));
if ~isempty(files)
    % Unzip files if they exist
    for i = 1:numel(files)
        zip_file = fullfile(data_dir, ['sub_' subject], files(i).name); % Corrected variable name
        fprintf('Image has not been unzipped; unzipping %s\n', zip_file);
        gunzip(zip_file);
    end
else
    fprintf('No images to unzip for subject %s\n', subject);
end

%%%%%%%%%%

% Specify scans and volumes
scans = spm_select('ExtFPList', fullfile(data_dir) ,'.*.nii$', 1:476);

clear matlabbatch

% Apply a 5 mm FWHM Gaussian smoothing kernel
matlabbatch{1}.spm.spatial.smooth.data = cellstr(scans);
matlabbatch{1}.spm.spatial.smooth.fwhm = [5 5 5];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

% Run job
spm_jobman('run', matlabbatch);

end