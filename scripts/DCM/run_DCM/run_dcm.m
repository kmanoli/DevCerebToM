function dcm

% FIRST-LEVEL DYNAMIC CAUSAL MODELLING
%
% SPM12
%
% This script was originally developed by Frank Van Overwalle's group and 
% edited by Katerina Manoli for this project. This analysis follows Van
% Overwalle, Van de Steen & Mariën's (2019; CABN) method.
%
% You should first specify models using the SPM GUI (DCM > specify) for one  
% specific participant (e.g., pp01). This initial DCM was specified for the
% first participant in each sample via the SPM12 GUI in the following way:

% - We included both conditions (ToM and pain)
% - We included the VOI timings provided by default in the GUI
% - We included the TE specified in the acquisition parameters in each
% study

% Model options window: 
% - Modulatory effects: bilinear
% - States per region: one
% - States per region: one
% - Stochastic effects: no
% - Centre input: no
% - Fit timeseries or CSD: timeseries

% Specification of connections windows:
% We wanted to specify all possible fixed connections between ROIs in 
% the first pop-up window. All these values are automatically turned on 
% for estimation of this (full) model - no need to select anything. It is
% important to also leave the following two pop-up windows (modulatory 
% connections)blank.

% This creates a *.mat file in the results subfolder of the participant 
% (e.g., the folder where the results of a first-level contrast analysis are). 
% Make a novel 'DCM' subfolder in your study folder (i.e., the folder where 
% all the participant's folders are located) and copy this *.mat file in it. 
% This script automatically uses these data of the first participant for 
% all the participants in the folder.
%
% This script looks for all models 'DCM_*.mat' in the DCM folder and
% estimates them. You can narrow down the number of models tested by
% specifying a model prefix that is inserted before the * in 'DCM_*.mat'
%
% After you have followed these steps, you can press 'Run'. Note that
% the model option used in the input box for this study was 1 (i.e., 
% experimental conditions ToM & pain placed first in the matrix).

% NOTE: 
% If you want to run more than 8 ROIs (e.g., 10), you have to change 
% the code in file "spm_dcm_specify_ui_m" in folder "spm12":
% In line 46, replace "8" by "10" (or a larger number if more ROIs are 
% required): 

% [P, sts] = spm_select([1 8],'^VOI.*\.mat$',{'select VOIs'},'',swd);

% so that it looks like this:
% [P, sts] = spm_select([1 10],'^VOI.*\.mat$',{'select VOIs'},'',swd);
%
%--------------------------------------------------------------------------

% create inputbox
answer = inputdlg({'Participant folder prefix (full folder for 1 ppt):', 'Start from participant:', 'results subfolder of each participant', 'Model name prefix (before * of DCM_*.mat)', 'Model option (experimental conditions first = 1 or last = 2; 0 = none)', 'Analysis (= time or csd)'}, 'Input', 1, ...
    {'PP', '1', 'Results_withAllTrials_New-art', 'Full7', '1', 'time'});
FolderPrefix = [answer{1}, '*']
StartingParticipant = str2num(answer{2})
RESsubfolder = answer{3}
ModelPrefix = [answer{4}, '*']
ModelOption = str2num(answer{5})
DCManalysis = answer{6}


%Make a list of all participant folders:
ppfolders = dir(FolderPrefix); 
pp = (size(ppfolders, 1))
    %empty ppdirectories
    ppdirectories = {};
    %for each participant
    for S1 = 1:pp
        try
            ppdirectories{S1}= cd(cd(ppfolders(S1).name));
        continue
        end
    end
    
%Make a list of all models:
mmfiles = dir(strcat('DCM/DCM_', ModelPrefix,'.mat'));
mm = (size(mmfiles, 1))
    %empty mmdirectories
    mmdirectories = {};
    %for each model
    for M1 = 1:mm
        try
            mmdirectories{M1,1} = mmfiles(M1).name;
        continue
        end
    end
    disp('Models (the 1st model is filled as full, see also BMR):')
    for M1 = 1:mm
        disp(['   ' mmdirectories{M1,1}])
    end
    
F = [];
for S2 = StartingParticipant:pp  
   
    dataDir = ppdirectories{S2}

    % change SPM.mat with the current images subfolder if needed
    ss = fullfile(dataDir, RESsubfolder, 'SPM.mat');
    load(ss, 'SPM');

    % IF YOUR PREPROCESSED FILES ARE NOT IN A PREPROCESSED SUBFOLDER, 
    % LEAVE "PREPROCESSED/" BLANK IN THE NEXT 10 LINES 
    imagefiles = dir([dataDir '/Preprocessed/swcur*.nii']);
    
    if strcmp(SPM.xY.VY(1,1).fname, [dataDir '/Preprocessed/' imagefiles(1).name]) == 0
        DirInfo = dir(ss);
        copyfile(ss, [ss '_' DirInfo.date(1:11)], 'f');
        disp(['Original SPM.mat with prior image paths copied into  ' ss '_' DirInfo.date(1:11)])
        imagefileschar = [];
        for i = 1: numel(imagefiles)
            SPM.xY.VY(i,1).fname = [dataDir '/Preprocessed/' imagefiles(i).name];
            imagefileschar = [imagefileschar; [dataDir '/Preprocessed/' imagefiles(i).name]];
        end 
        SPM.xY.P = imagefileschar;
        disp(['Adjusted SPM.mat with current image paths created as ' ss])
        save(ss);
    end
    
    
    % Load first model
    load(fullfile(pwd, 'DCM', mmdirectories{1}));

    % Load ROIs and check whether they exist for this participant
    ROInum = DCM.n; 
    for i = 1: ROInum
        ROIlist{i} = DCM.xY(1, i).name;
        ROIfile = strcat(dataDir, '/', RESsubfolder, '/VOI_', ROIlist{i}, '_1.mat');
        ROIfiles{i, 1} = ROIfile;
        NotAllROIs = 0;
        if isempty(dir(ROIfile))
            NotAllROIs = 1;
        end
    end 
    if (NotAllROIs == 1); 
        disp('This participant does not have all required VOIs')
        break; 
    end

    % Load all DCM models 
    % Load experimental data for selected conditions
    for jj = 1: mm
        load(fullfile(pwd, 'DCM', mmdirectories{jj}));      
        % load ROI time series 
        DCM.Y.y = [];
        for rr = 1: ROInum
            ROIfile = ROIfiles{rr, 1};
            load(ROIfile,'xY');                     % xY = each ROI
            DCM.xY(rr) = xY;                        %
            DCM.Y.y(:,rr)  = DCM.xY(rr).u;          % xY.u = time series of each ROI
        end
        % load scans with onsets for each condition that is included
        DCM.U.u = [];
        CONDnum = 0;
        for u1 = 1:length(SPM.Sess.U)
            for u2 = 1:length(DCM.U.name)
                try
                    if (SPM.Sess.U(u1).name{1} == DCM.U.name{u2})
                    DCM.U.u = [DCM.U.u SPM.Sess.U(1,u1).u(33:end,1)];  
                    CONDnum  = CONDnum + 1;  
                    end
                end
            end
        end
        % load to set the correct size of arrays
        DCM.v    = length(DCM.xY(1).u);             % number of time points
        DCM.Y.Q  = spm_Ce(ones(1,DCM.n)*DCM.v);
        DCM.Y.X0 = DCM.xY(1).X0;                    
        % Make first model a full model
        if jj == 1
            for i = 1:ROInum                        % Number of VOIs
                for j = 1:ROInum
                    for cond = 1:CONDnum            % Number of conditions
                        DCM.a(i,j) = 1;
                        DCM.b(i,j,cond) = 1;
                        DCM.c(i,cond) = 1;
                    end
                end
            end
        end
        
        if isempty(dir(fullfile(dataDir, 'DCM')));
            mkdir(fullfile(dataDir, 'DCM'));
        end
        save(fullfile(dataDir, 'DCM', mmdirectories{jj}),'DCM');
    end
    
    % DCM restructuring into less redundant model (option: experimental conditions first = 1, last = 2)
    if ModelOption > 0
        % specify "all" control values
        Uall = DCM.U.u(:,1);
        Tall = DCM.U.name{1}(1:3);
        Call = DCM.c(:,1);
        for ii = 2:CONDnum
            Uall = Uall + DCM.U.u(:,ii);
            Tall = [Tall '+' DCM.U.name{ii}(1:3)];
            Call = max(Call, DCM.c(:,ii));
        end
        % structuring DCM 
        for ii = 1:CONDnum - 1
        % shift experimental values one step down if they come last  
        if ModelOption == 2
            DCM.U.u(:,ii)  = DCM.U.u(:,ii + 1);
            DCM.U.name{ii} = DCM.U.name{ii + 1};
            DCM.b(:,:,ii)  = DCM.b(:,:,ii + 1);        
        end
        DCM.c(:,ii)    = 0;                % set input of experimental conditions to zero
        end
        % set "all" control values last
        DCM.U.u(:,CONDnum)  = Uall; 
        DCM.U.name{CONDnum} = Tall;
        DCM.b(:,:,CONDnum)  = 0;           % set modulation of "all" conditions to zero
        DCM.c(:,CONDnum)    = Call;        % set input of "all" conditions to maximum (of 1 and 0)   
        save(fullfile(dataDir, 'DCM', mmdirectories{jj}),'DCM');
    end
    
    % Model name
    ModelName = mmdirectories{1};
    [p, n, e] = fileparts(ModelName);
    n = n(5:length(n));

    % DCM Estimation
    clear matlabbatch
    %matlabbatch{1}.spm.dcm.spec.fmri.regions.dcmmat = cellstr(construct_list([dataDir '/DCM'], mmdirectories, '', 0));
    %matlabbatch{1}.spm.dcm.spec.fmri.regions.voimat = ROIfiles;
    %matlabbatch{2}.spm.dcm.spec.fmri.inputs.dcmmat = cellstr(construct_list([dataDir '/DCM'], mmdirectories, '', 0));
    %matlabbatch{2}.spm.dcm.spec.fmri.inputs.spmmat = {[dataDir '/' RESsubfolder '/SPM.mat']};
    %matlabbatch{2}.spm.dcm.spec.fmri.inputs.session = 1;
    %matlabbatch{2}.spm.dcm.spec.fmri.inputs.val = {1};    %??
    matlabbatch{1}.spm.dcm.estimate.dcms.model.dcmmat = cellstr(construct_list([dataDir '/DCM'], mmdirectories, '', 0));
    matlabbatch{1}.spm.dcm.estimate.output.single.dir = {[dataDir, '/DCM']};
    matlabbatch{1}.spm.dcm.estimate.output.single.name = [n '+nested'];
    matlabbatch{1}.spm.dcm.estimate.est_type = 1;
    matlabbatch{1}.spm.dcm.estimate.fmri.analysis = DCManalysis;
    
    batchfilename = [ppdirectories{S2} '/DCM/' 'batch_subject' num2str(S2) '_dcm.mat']; 
    save(batchfilename, 'matlabbatch');
    
    try % in case there is no convergence
    spm_jobman('run', matlabbatch);
    
%--------------------------------------------------------------------------
% Collate subject GCMs in group GCM structure

    GCMsubject = load(fullfile(dataDir, 'DCM', ['GCM_' n '+nested.mat']),'GCM');
    disp(strcat('Writing subject ', num2str(S2), ' DCM in a group GCM... ', pwd, '/DCM/GCM_', n, '+nested.mat'));
        % Load models
        for M2 = 1:mm
            GCMsample{S2, M2} = GCMsubject.GCM{1, M2};
        end
        GCM = GCMsample;
    % Write group GCM structure
    save(fullfile(pwd, 'DCM', ['GCM_' n '+nested.mat']),'GCM', spm_get_defaults('mat.format'));

%--------------------------------------------------------------------------
% Write out in Excel & F.mat

% Load models
for M2 = 1:mm
    ModelName = mmdirectories{M2};

    % Load data
    %DCM = (GCM{1,M2});
    clear DCM
    load(fullfile(dataDir, 'DCM', [ModelName]))
    % For testing excel file
    %load(fullfile(dataDir, 'DCM', ['Copy_of_' ModelName]))
    
    ROInum = DCM.n; 
    ROIlist=[];
    for i = 1: ROInum
        ROIlist{i} = DCM.xY(1, i).name;
    end
    [p, n, e] = fileparts(ModelName);
    xlsfilename = fullfile(pwd, strcat('DCM/', n , '.xls'));
    disp(strcat('Writing results on excel file... ',xlsfilename));
    [p, n, e] = fileparts(ppdirectories{S2});
    
    % Define participant's column by alphabetical Excel reference
    xlChr = '';
    xlCol = (S2 - 1) * max(ROInum, CONDnum) + 3;
    for i = 1:(max(ROInum, CONDnum) * pp);
        if (xlCol > 26);
            xlChr = char(i + 64);
            xlCol = xlCol - 26;
        end
        cond = char(xlCol + 64);
    end
    xlChr = strcat(xlChr, cond); 
    xlswrite(xlsfilename, {n}, 1, strcat(xlChr, '1'));
    xlswrite(xlsfilename, ROIlist, 1, strcat(xlChr, '2'));            

    % Write matrices
    lbl = 'A: Fixed - estimates';
    disp(lbl);
        xlRow = 0 * (ROInum + 1) + 2;
        if (S2 == 1); xlswrite(xlsfilename, {lbl}, 1, strcat('A',num2str(xlRow + 1))); end
        if (S2 == 1); xlswrite(xlsfilename, ROIlist', 1, strcat('B',num2str(xlRow + 2))); end
        xlswrite(xlsfilename, DCM.Ep.A, 1, strcat(xlChr, num2str(xlRow + 2)));
    lbl = 'A: Fixed - probabilities';
    disp(lbl);
        xlRow = 1 * (ROInum + 1) + 2;
        if (S2 == 1); xlswrite(xlsfilename, {lbl}, 1, strcat('A',num2str(xlRow + 1))); end
        if (S2 == 1); xlswrite(xlsfilename, ROIlist', 1, strcat('B',num2str(xlRow + 2))); end
        xlswrite(xlsfilename, DCM.Pp.A, 1, strcat(xlChr, num2str(xlRow + 2)));
 
    try % in case of resting state models

    lbl = 'B: Modulatory - estimates';
    disp(lbl);
    for i = 1:CONDnum;
        cond = DCM.U.name{i};
        xlRow = (2 + i - 1) * (ROInum + 1) + 3;
        if (S2 == 1); xlswrite(xlsfilename, {lbl}, 1, strcat('A',num2str(xlRow + 1))); end
        if (S2 == 1); xlswrite(xlsfilename, {cond}, 1, strcat('A',num2str(xlRow + 2))); end
        if (S2 == 1); xlswrite(xlsfilename, ROIlist', 1, strcat('B',num2str(xlRow + 2))); end
        xlswrite(xlsfilename, DCM.Ep.B(:,:,i), 1, strcat(xlChr, num2str(xlRow + 2)));
    end
    lbl = 'B: Modulatory - probabilities';
    disp(lbl);
    for i = 1:CONDnum;
        cond = DCM.U.name{i};
        xlRow = (CONDnum + 2 + i - 1) * (ROInum + 1) + 3;
        if (S2 == 1); xlswrite(xlsfilename, {lbl}, 1, strcat('A',num2str(xlRow + 1))); end
        if (S2 == 1); xlswrite(xlsfilename, {cond}, 1, strcat('A',num2str(xlRow + 2))); end
        if (S2 == 1); xlswrite(xlsfilename, ROIlist', 1, strcat('B',num2str(xlRow + 2))); end
        xlswrite(xlsfilename, DCM.Pp.B(:,:,i), 1, strcat(xlChr, num2str(xlRow + 2)));
    end
    
    lbl = 'C: Direct input - estimates';
    disp(lbl);
    cond = []; for i = 1:CONDnum; cond = cellstr([cond; DCM.U.name{i}]); end
        xlRow = ((2 * CONDnum) + 2) * (ROInum + 1) + 4;
        xlswrite(xlsfilename, cond', 1, strcat(xlChr, num2str(xlRow + 1)));
        if (S2 == 1); xlswrite(xlsfilename, {lbl}, 1, strcat('A',num2str(xlRow + 1))); end
        if (S2 == 1); xlswrite(xlsfilename, ROIlist', 1, strcat('B',num2str(xlRow + 2))); end
        xlswrite(xlsfilename, DCM.Ep.C, 1, strcat(xlChr, num2str(xlRow + 2)));
    lbl = 'C: Direct input - probabilities';
    disp(lbl);
        xlRow = ((2 * CONDnum) + 3) * (ROInum + 1) + 4;
        if (S2 == 1); xlswrite(xlsfilename, {lbl}, 1, strcat('A',num2str(xlRow + 1))); end
        if (S2 == 1); xlswrite(xlsfilename, ROIlist', 1, strcat('B',num2str(xlRow + 2))); end
        xlswrite(xlsfilename, DCM.Pp.C, 1, strcat(xlChr, num2str(xlRow + 2)));

    catch
        % this is a resting state model
    end
        
    lbl = 'F: log-evidence';
    disp(lbl);
        xlRow = ((2 * CONDnum) + 4) * (ROInum + 1) + 5;
        if (S2 == 1); xlswrite(xlsfilename, {lbl}, 1, strcat('A',num2str(xlRow + 1))); end
        xlswrite(xlsfilename, DCM.F, 1, strcat(xlChr, num2str(xlRow + 2)));
        
    %Store log-evidence in F structure
    %disp(strcat('Writing log-evidence f on mat-file... ', pwd, '/DCM/F.mat'));
    %F(S2 - StartingParticipant + 1, M2) = DCM.F;
    %save(fullfile(pwd,['DCM/F.mat']),'F', spm_get_defaults('mat.format'));
        
    end
    
    end % end of try
    
end

% rename with random postscript and model name to avoid overwriting
rng('shuffle');
r = num2str(randi(1000, 1));
    ModelName = mmdirectories{1};
    [p, n, e] = fileparts(ModelName);
    n = n(5:length(n));
movefile(fullfile(pwd, 'DCM', ['GCM_' n '+nested.mat']), ...
    fullfile(pwd, 'DCM', ['GCM_' n '+nested_' r '.mat']), 'f')
movefile(fullfile(pwd,['DCM/F.mat']), ...
    fullfile(fullfile(pwd,['DCM/F_' r '.mat'])), 'f')
for M2 = 1:mm
    ModelName = mmdirectories{M2};
    [p, n, e] = fileparts(ModelName);
    xlsfilename = fullfile(pwd, strcat('DCM/', n , '.xls'));
    movefile(xlsfilename, ...
        [xlsfilename(1:length(xlsfilename) - 4) '_' r '.xls'], 'f')
end

end
    
function files = construct_list(dir, filelist, prefix, listtype)
%if listtype==1
%    countlist = num2str(ones(size(filelist, 1), 1));
%else
%    countlist = num2str((1:size(filelist, 1))');
%end
%files = strcat(dir, filesep, prefix, filelist, ',', countlist);
files = strcat(dir, filesep, prefix, filelist);
end
 
 

