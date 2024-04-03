function runPEB

% GROUP-LEVEL DYNAMIC CAUSAL MODELLING
%
% SPM12
%
% This script was originally developed by Frank Van Overwalle's group and 
% edited by Katerina Manoli for this project. This analysis follows Van
% Overwalle, Van de Steen & Mariën's (2019; CABN) method.
%
% This script runs a group-level Parametric Empirical Bayes (PEB) model, 
% which makes it possible to evaluate group effects and between-subjects 
% variability on connectivity parameters (Friston et al., 2016).
%
% In order to test differences between groups, you can provide an .xls
% file with covariates. The first row should denote the variable name 
% (e.g., ToM task performance), and the next cells determine the contrasts 
% for each row (i.e., participant). Here, we coded ToM passers as 1 and 
% non-passers as -1. 
%
% After you have followed these steps, you can press 'Run'. Note that
% the model option used in the input box for this study was 1 (i.e., 
% experimental conditions ToM & pain placed first in the matrix).
%
% NOTE:
% The script saves the following files:
% - A batch file with prefix PEBsearch (e.g., PEBsearch*.mat),which you can 
% inspect in case of problems.
% - An Excel file with prefix BMA_PEB (e.g., BMA_PEB*.xls), with the posterior
% estimates (Ep) and posterior probabilities (Pp) in an easy-to-read format 
% for all models and all covariates/contrasts. Matrices A, B  denote  
% fixed and modulatory connections, respectively. In case covariates were
% included, the *second* A, B matrices denote the differences between groups
% for the specified contrast(s).
%--------------------------------------------------------------------------

clc;clear;

% create inputbox
answer = inputdlg({'Number of DCM models', 'Number of Covariates (0 = none)', 'group GCM*.mat name', 'Covariates*.xls name', 'Number of Conditions (including control)', 'Model Option (experimental conditions first = 1 or last = 2; 0 = none)', 'Analysis (= time or csd)'}, 'Input', 1, ...
    {'1', '0', 'GCM_Full7+nested_358.mat', ' ', '2', '1', 'time'});
%    {'1', '1', 'GCMs_full+nested.mat', 'Covariates.xls', '5', '2', 'time'});
%    {'1', '2', 'GCM_groups_full+nested_Patients+Healthy_DN.mat', 'Covariates_ALL.xls', '0', '0', 'csd'});
%    {'1', '2', 'GCM_groups_full+nested_ALL_DN.mat', 'Covariates_ALL.xls', '0', '0', 'csd'});
%    {'1', '3', 'GCM_full+nested.mat', 'Covariates_n=73.xls', '2', '2', 'time'});
%    {'1', '1', 'GCM_full+nested.mat', 'Covariates_SocialConsistent.xls', '2', '2', 'time'});
%    {'1', '4', 'GCM_full+nested.mat', 'Covariates_perStudy.xls', '2', '2', 'time'});
%    {'1', '0', 'GCM_full+nested.mat', ' ', '2', 'time'});

mm          = str2num(answer{1})
COVnum      = str2num(answer{2})
GCMfilename = answer{3};
COVfilename = answer{4};
CONDnum     = str2num(answer{5})
ModelOption = str2num(answer{6})
DCManalysis = answer{7}


cd ('DCM')

GCM = load(fullfile(pwd, GCMfilename),'GCM')%;
DCMfilename = char(GCM.GCM(1,1))%;
DCM = load(fullfile(DCMfilename),'DCM')%;

COVsubj = size(GCM.GCM,1)%;

if COVnum > 0

    % read labels on top of 'Covariates' file (first row)
    [num, txt, head] = xlsread(COVfilename, 'A1:Z1000');
    [nrow, ncol] = size(num);
    if nrow ~= COVsubj; 
        disp(['ERROR: Number of participants (n=' num2str(COVsubj) ') do not match with Covariance Design file (n=' num2str(nrow) ')']); 
        disp(['ERROR: Avoid all numbers in the first row of the Covariance Design excel file!']); 
        return
    end
    xlChr = strcat(char(COVnum + 64));
    [num, cov_names] =  xlsread(COVfilename, ['A1:' xlChr, '1']);
    cov_design = xlsread(COVfilename, ['A2:' xlChr, num2str(COVsubj + 1)]);
end 

%--------------------------------------------------------------------------
% Run the batch 

% !! only tested on 1 model

if COVnum > 0
    filename = [GCMfilename(5:length(GCMfilename) - 4) '_' COVfilename(1:length(COVfilename) - 4)];
else
    filename = [GCMfilename(5:length(GCMfilename) - 4) '_NoCovariates']%;
end

for M2 = 1:mm
    matlabbatch{1}.spm.dcm.peb.specify.name = filename%;
    matlabbatch{1}.spm.dcm.peb.specify.model_space_mat = cellstr(construct_list(pwd, GCMfilename, '', 0))%;
    matlabbatch{1}.spm.dcm.peb.specify.dcm.index = 1;
    if COVnum > 0 
        matlabbatch{1}.spm.dcm.peb.specify.cov.design_mtx.cov_design = cov_design;
        matlabbatch{1}.spm.dcm.peb.specify.cov.design_mtx.name = cov_names;
    end;
    if ~strcmpi(DCManalysis, 'csd')    
    matlabbatch{1}.spm.dcm.peb.specify.fields.custom = {
                                                    'A'
                                                    'B'
                                                    'C'
                                                    }'; 
    end;
    matlabbatch{1}.spm.dcm.peb.specify.priors_between.ratio = 16;
    matlabbatch{1}.spm.dcm.peb.specify.priors_between.expectation = 0;
    matlabbatch{1}.spm.dcm.peb.specify.priors_between.var = 0.0625;
    matlabbatch{1}.spm.dcm.peb.specify.priors_glm.group_ratio = 1;
    matlabbatch{1}.spm.dcm.peb.specify.show_review = 0;
    matlabbatch{2}.spm.dcm.peb.reduce_all.peb_mat(1) = cfg_dep('Specify / Estimate PEB: PEB mat File(s)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','peb_mat'))%;
    matlabbatch{2}.spm.dcm.peb.reduce_all.model_space_mat = cellstr(construct_list(pwd, GCMfilename, '', 0))%;
    matlabbatch{2}.spm.dcm.peb.reduce_all.nullpcov = 0.0625;
    matlabbatch{2}.spm.dcm.peb.reduce_all.show_review = 1;

    PEBfilename = ['PEBsearch_' filename]%;
    batchfilename = [PEBfilename '.mat']%; 
    %%%save(batchfilename, 'matlabbatch');
    spm_jobman('run', matlabbatch);
    
%--------------------------------------------------------------------------
% Write out in Excel 

GCM = load(fullfile(pwd, GCMfilename),'GCM')%;
DCMfilename = char(GCM.GCM(1,1))%;
DCM = load(fullfile(DCMfilename),'DCM')%;

%BMAfilename = ['BMA_PEB_' filename '.mat'];
BMAfilename = ['BMA_search_PEB_' filename '.mat'];
BMA = load(fullfile(pwd, BMAfilename),'BMA');

%xlsfilename = fullfile(pwd, ['BMA_PEB_' filename '.xls']);
xlsfilename = fullfile(pwd, ['BMA_search_PEB_' filename '.xls']);
disp(strcat('Writing results on excel file... ',xlsfilename));

% !! only tested on 1 model
% !! only tested on ModelOption > 0

% Load models
for M2 = 1:mm 
    n = ['Model_' num2str(M2)];

    % Load data
    %ROInum = size(GCM.GCM{1,1}.xY, 2); 
    ROInum = size(DCM.DCM.xY, 2); 
    ROIlist=[];
    for i = 1:ROInum
        %ROIlist{i} = GCM.GCM{1,1}.xY(1,i).name;
        ROIlist{i} = DCM.DCM.xY(1,i).name;
    end
    
    COVnum = size(BMA.BMA.Xnames, 2); 
    COVlist=[];
    for i = 1:COVnum
        COVlist{i} = BMA.BMA.Xnames{1,i};
    end    
    
    % Define model's column by alphabetical Excel reference
    xlChr = '';
    xlCol = (M2 - 1) * max(ROInum, CONDnum) + 3;
    for i = 1:(max(ROInum, CONDnum) * mm);
        if (xlCol > 26);
            xlChr = char(i + 64);
            xlCol = xlCol - 26;
        end
        cond = char(xlCol + 64);
    end
    xlChr = strcat(xlChr, cond); 
            for r = 1:ROInum 
                xlCh{r} = char(xlChr + r - 1);
            end;    

    % Write Group mean and Covariates
    %xlswrite(xlsfilename, {n}, 1, strcat(xlChr, '1'));            
    %xlswrite(xlsfilename, ROIlist, 1, strcat(xlChr, '2'));
    
    writecell({n},xlsfilename, 'Range', strcat(xlChr, '1'));
    writecell(ROIlist,xlsfilename, 'Range', strcat(xlChr, '2'));

    for covar = 1:COVnum
        disp (' ');
        disp([char(n), ' - ', char(COVlist(covar))]);

        % Write matrices
        % time Analysis
        Anum = 1;
        Bnum = CONDnum;
        Cnum = CONDnum;
        if ModelOption > 0; 
            Bnum = CONDnum - 1;
            Cnum = 1;
        end
        nn = Anum + Bnum + Cnum;
        % csd Analysis
        if strcmpi(DCManalysis, 'csd'); nn = Anum; Cnum = 0; end    

        % Writing matrices
        for matrix = 1:nn
            if matrix == 1; lbl = 'A: Fixed - estimates'; end
            if matrix  > 1; lbl = 'B: Modulatory - estimates'; end
            if matrix  > Anum + Bnum; lbl = 'C: Input - estimates'; end
            if matrix == 1; cond = ''; end
            %if matrix  > 1; cond = char(GCM.GCM{1, 1}.U.name{matrix - 1});
            if matrix  > 1; cond = char(DCM.DCM.U.name{matrix - 1});
            end
            if matrix  > Anum + Bnum
                rr = (covar - 1)*((Anum + Bnum)*(ROInum^2) + Cnum*ROInum) + (Anum + Bnum)*(ROInum^2) + (matrix - Anum - Bnum - 1)*ROInum; 
            else
                rr = (covar - 1)*((Anum + Bnum)*(ROInum^2) + Cnum*ROInum) + (matrix - 1)*(ROInum^2);
            end   
            if strcmpi(DCManalysis, 'csd'); 
                rr = (covar - 1)*((ROInum^2)*nn) + (matrix - 1)*(ROInum^2);
            end    
            disp(lbl);
                xlRow = (covar - 1)*((Anum + Bnum + Cnum)*(ROInum + 1)*2) + (matrix - 1)*(ROInum + 1)*2 + 2; 
                if (M2 == 1)
                    %xlswrite(xlsfilename, [{n} COVlist(covar)], 1, strcat('D',num2str(xlRow + 1))); 
                    %xlswrite(xlsfilename, {'from'}, 1, 'B1'); 
                    %xlswrite(xlsfilename, {'to'}, 1, 'A2'); 
                    %xlswrite(xlsfilename, {lbl}, 1, strcat('A',num2str(xlRow + 1))); 
                    %xlswrite(xlsfilename, {cond}, 1, strcat('A',num2str(xlRow + 2))); 
                    %xlswrite(xlsfilename, ROIlist', 1, strcat('B',num2str(xlRow + 2))); 
                    
                    writecell([{n} COVlist(covar)],xlsfilename,'Range', strcat('D',num2str(xlRow + 1)));
                    writecell({'from'},xlsfilename,'Range','B1');
                    writecell({'to'},xlsfilename,'Range','A2');
                    writecell({lbl},xlsfilename,'Range', strcat('A',num2str(xlRow + 1)));
                    writecell({cond},xlsfilename,'Range', strcat('A',num2str(xlRow + 2)));
                    writecell(ROIlist',xlsfilename,'Range', strcat('B',num2str(xlRow + 2)));
                end
                if lbl(1:1) ~= 'C'; num = ROInum; else num = Cnum; end
                % !! Only tested for Model Option = 2; !!! check for ModelOption = 0
                for r = 1:num
                    %xlswrite(xlsfilename, BMA.BMA.Ep(rr + (r - 1)*ROInum + 1:rr + r*ROInum), 1, strcat(char(xlCh{r}), num2str(xlRow + 2)));
                    writematrix(BMA.BMA.Ep(rr + (r - 1)*ROInum + 1:rr + r*ROInum),xlsfilename,'Range',strcat(char(xlCh{r}), num2str(xlRow + 2)));
                end
            if matrix == 1; lbl = 'A: Fixed - probabilities'; end
            if matrix  > 1; lbl = 'B: Modulatory - probabilities'; end
            if matrix  > Anum + Bnum; lbl = 'C: Input - Probabilities'; end
            disp(lbl);
                xlRow = (covar - 1)*((Anum + Bnum + Cnum)*(ROInum + 1)*2) + (matrix - 1)*(ROInum + 1)*2 + ROInum + 3; 
                if (M2 == 1); 
                    %xlswrite(xlsfilename, {lbl}, 1, strcat('A',num2str(xlRow + 1))); 
                    %xlswrite(xlsfilename, {cond}, 1, strcat('A',num2str(xlRow + 2))); 
                    %xlswrite(xlsfilename, ROIlist', 1, strcat('B',num2str(xlRow + 2))); 
                    writecell({lbl},xlsfilename,'Range', strcat('A',num2str(xlRow + 1)));
                    writecell({cond},xlsfilename,'Range', strcat('A',num2str(xlRow + 2)));
                    writecell(ROIlist',xlsfilename,'Range', strcat('B',num2str(xlRow + 2)));
                end
                if lbl(1:1) ~= 'C'; num = ROInum; else num = Cnum; end
                for r = 1:num
                    %xlswrite(xlsfilename, BMA.BMA.Pp(rr + (r - 1)*ROInum + 1:rr + r*ROInum), 1, strcat(char(xlCh{r}), num2str(xlRow + 2)));
                    writematrix(BMA.BMA.Pp(rr + (r - 1)*ROInum + 1:rr + r*ROInum),xlsfilename,'Range',strcat(char(xlCh{r}), num2str(xlRow + 2)));

                end
        end
    end
 
end

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
    
 

