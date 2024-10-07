%% Script to run gPPI analysis on Responsibility study data.
% As we use spm8 to run the PPPI and spm12 to estimate the contrasts, it's
% recommended to exit Matlab between model estimation and contrast calculation.
%
% J Schultz 2019-2024

%% Preparation
clear, close all
addpath(genpath('/Applications/PPPIv13.1'))

%% Flags
deletePreviousPPI = 0;
runPPI = 0;
runContrasts = 1;
deletePreviousContrasts = 1;
displayContrastCheck = 0;

%% Settings
codeFolder = '/Users/johannesschultz/sciebo/CURRENT WORK/Responsibility/Code_minimal/'; % just used for error reporting
mainFolder = '/Volumes/LaCie_RAID_A/Documents/exp/SODEC/Responsibility/FFX/'; % in which the FFXname folder is
FFXname = 'main/'; % the model used for all the main results in the paper. This folder contains the individual FFX folders

% Conditions to be investigated, from the FFX model !! need to be identical 
% to condition names in FFX, but no mention of session or basis function.
FFXconds = {'0'...
  'Decision safe active'...
  'Decision safe passive'...
  'Decision safe n-soc'...
  'Decision risky active'...
  'Decision risky passive'...
  'Decision risky n-soc'...
  } ;

% Participants
persons_to_do = [ 44 55 62 63 65 69 70 71 74 77    79 81 83 84 85    90 92 93 94 95 96 98 101 102 103 104 106 110 120     124 125 126 129 136 155 166 171 173 174 175 ]; % all participants except those without data in OFC (87 and 176), people with only 1 session (121) or aborted Sess 1 (78), too few "risky" choices (87), missing risky outcome other pos social and self neg solo in Sess 1 (177). Still in are 62 (only 1 trial of outcome other neg partner in Sess 1), 65 and 120, missing outcome self neg solo in both sess (65) and Sess 1 (120)
% persons_to_do = [ 44 ]; % test participant
persons_to_do = [ 44 55 62 63 65 69 70 71 74 77    79 81 83 84 85    90 92 93 94 95 96 98 101 102 103 104 106 110 120     124 125 126 129 136 155 166 171 173 174 175 ]; % all participants except those without data in OFC (87 and 176), people with only 1 session (121) or aborted Sess 1 (78), too few "risky" choices (87), missing risky outcome other pos social and self neg solo in Sess 1 (177). Still in are 62 (only 1 trial of outcome other neg partner in Sess 1), 65 and 120, missing outcome self neg solo in both sess (65) and Sess 1 (120)

% Define the seed ROIs
% The way I ran the original analysis: get ROIs by listing: ROIs = dir([ROIfileDir '*.nii']); % I use functional ROIs with the same resolution and mask
%   ROIs{1}     = [mainFolder 'PPIseeds/GuiltEffect_InsulaLeftMaskedByRisky>SafeOutcome_0p001u.nii']; % The new, small left anterior insula ROI
%   PPInames{1} = 'PPI_GuiltEffect_InsulaLeft';
%   % ROIs{1} = [mainFolder 'PPIseeds/guiltEffect_0p05FWE_SVC_aIns.nii']; % This is the latest version of the left aIns cluster. Use this??

% ROIs{1}     = [mainFolder 'PPIseeds/STSleft_otherRPE_actPass_0p05FWEclust.nii']; % The STS cluster identified in the comp model-based analysis
% PPInames{1} = 'PPI_STSleft_sRPEsocPartner';

ROIs{1}     = [mainFolder 'PPIseeds/insulaLeft_outcomeRiskySafe_0p001FWEvox.nii']; % The insula cluster identified in outcomeRiskySafe
ROInames{1} = 'insLeft_outcRiskySafe'; % the PPI model will named "PPI_" followed by this ROIname, and located within the FFX it is based on

%% Run PPI FFX models
if runPPI
  spm8
  addpath(genpath('PPPIv13.1'))
  
  for person_i = 1 : length( persons_to_do )
    personNo = persons_to_do(person_i);
    personName = sprintf('SODEC_FMRI_%02d',personNo);
    FFXdir = [mainFolder FFXname personName];
    if deletePreviousPPI
      owd = pwd;
      cd(FFXdir)
      try rmdir('PPI*','s'); end % in case they already exist, PPPI will create models in these directories but not estimate them!
      try delete gppi* ; end
      cd(owd)
    end
    
    for r = 1%:length(ROIs)
      
      clear P
      
      P.subject = ['gppi_' personName '_Seed_' ROInames{r}];
      P.directory = FFXdir;
      
      P.VOI = ROIs{r}; % [ROIfileDir ROIs(r).name];
      
      P.Region = ROInames{r};
      
      P.Estimate = 1;
      P.contrast = 0;
      P.extract = 'eig';
      
      P.Tasks = FFXconds; % get the correct condition names!
      
      P.Weights = [];
      P.analysis = 'psy';
      P.method = 'cond';
      P.CompContrasts = 0;
      P.Weighted = 0;
      P.ConcatR = 0;
      P.preservevarcorr = 0;
      P.wb = 0;
      
      PPPI(P)
      % to get notified of the progress:
      logfiles = dir([P.directory '/*.log']); logfile = logfiles(end);
      copyfile([logfile.folder filesep logfile.name],[codeFolder '/fMRI_gPPI_logs/' logfile.name]);
    end
    fclose('all'); % to avoid "Too many files open" error
  end % persons
end % runPPI
if runPPI; disp('Ran all models, do the contrasts now - see later in the script'); end

%% Run contrasts
% As we use spm8 to run the PPPI and spm12 to estimate the contrasts, it's
% recommended to exit Matlab between model estimation and contrast
% calculation.

if runContrasts
  disp('Calculating contrasts...')
  if isempty(which('spm')) || ~strcmpi(spm('ver'),'SPM12'); spm12; spm( 'defaults', 'fmri' ); end
  spm_jobman( 'initcfg' );
  % persons_to_do = [ 44 55 62 63 65 69 70 71 74 77    79 81 83 84 85    90 92 93 94 95 96 98 101 102 103 104 106 110 120     124 125 126 129 136 155 166 171 173 174 175 ]; % all participants except those without data in OFC (87 and 176), people with only 1 session (121) or aborted Sess 1 (78), too few "risky" choices (87), missing risky outcome other pos social and self neg solo in Sess 1 (177). Still in are 62 (only 1 trial of outcome other neg partner in Sess 1), 65 and 120, missing outcome self neg solo in both sess (65) and Sess 1 (120)
  % persons_to_do = [ 44 ]; % test
  persons_onlyScanned1 = [ 59 64 ]; % have only 1 sess (need separate con estimation)
  persons_onlyUseSess1 = [  ]; % has 2 sessions but only sess 1 is useable for contrasts
  persons_onlyUseSess2 = [  ]; % has 2 sessions but only sess 1 is useable for contrasts
  
  % 1-session contrast values
  conIni = zeros(1,32);
  i=1; conNames{i} = 'Decision safe social';    conVecs{i} = conIni; conVecs{i}(20) = 1;
  i=2; conNames{i} = 'Decision safe partner';   conVecs{i} = conIni; conVecs{i}(21) = 1;
  i=3; conNames{i} = 'Decision safe solo';      conVecs{i} = conIni; conVecs{i}(22) = 1;
  i=4; conNames{i} = 'Decision risky social';   conVecs{i} = conIni; conVecs{i}(23) = 1;
  i=5; conNames{i} = 'Decision risky partner';  conVecs{i} = conIni; conVecs{i}(24) = 1;
  i=6; conNames{i} = 'Decision risky solo';     conVecs{i} = conIni; conVecs{i}(25) = 1;
  i=7; conNames{i} = 'Risky>SafeXSocial>Solo';  con = conIni; con([22 23]) = 1; con([20 25]) = -1; conVecs{i} = con;
  i=8; conNames{i} = 'Risky>SafeXSolo>Social';  con = conIni; con([20 25]) = 1; con([22 23]) = -1; conVecs{i} = con;
  
  for person_i = 1 : length( persons_to_do )
    clear matlabbatch
    for r = 1:length(ROInames)
      PPIname = ['PPI_' ROInames{r}];
      personNo = persons_to_do(person_i);
      personName = sprintf('SODEC_FMRI_%02d',personNo);
      if ismember(personNo,persons_onlyUseSess1) % make a contrast weighting only session 1
        sessMult = [1 0];
      elseif ismember(personNo,persons_onlyUseSess2) % make a contrast weighting only session 2
        sessMult = [0 1];
      elseif ismember(personNo,persons_onlyScanned1) % make a contrast for only 1 session
        sessMult = [1];
      else, sessMult = [1 1];
      end
      disp([personName ,' PPI #' num2str(r) ': ' PPIname])
      FFXdir = [mainFolder FFXname personName filesep PPIname filesep];
      if exist(FFXdir,'dir')
        temp = load([FFXdir 'SPM.mat']); est = spm_SpUtil('IsCon',temp.SPM.xX.X);
      else
        disp(['Could not load SPM.mat from: ' FFXdir ]); return
      end
      try
        clear matlabbatch
        matlabbatch{1}.spm.stats.con.spmmat = { [FFXdir 'SPM.mat'] };
        for c = 1:length(conNames)
          % create contrast vector for all used sessions
          conVec = [kron(sessMult,conVecs{c}) zeros(1,length(sessMult))];
          conVecWorks = conVec.*est; % if sum(conVec.*est) == 0 % then at least the contrast is balanced, this was the case for subjects 98 and 175, so I adapted the contrast for them and used that
          if displayContrastCheck
            disp(conNames{c})
            disp('weighted positive:')
            char(temp.SPM.xX.name(conVecWorks>0))
            disp('weighted negative:')
            char(temp.SPM.xX.name(conVecWorks<0))
          end
          %                   if length(intersect(find(est),find(conVec))) == length(find(conVec)) % then all the betas needed are estimable
          %                     matlabbatch{1}.spm.stats.con.consess{cc+1}.tcon.name = conNames{c};
          %                     matlabbatch{1}.spm.stats.con.consess{cc+1}.tcon.convec = conVecWorks;
          %                     matlabbatch{1}.spm.stats.con.consess{cc+1}.tcon.sessrep = 'none';
          %                   cc = cc+1;
          matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = conNames{c};
          matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = conVecWorks;
          matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
          %                   else
          %                     fprintf('Contrast %d, %s not estimable\n',c,conNames{c})
          %                     figure; bar([conVec; est]'); legend({'Contrast','Estimable regressors'},'location','SouthWest')
          %                     keyboard
          %                   end
        end
        matlabbatch{1}.spm.stats.con.delete = deletePreviousContrasts;
        spm_jobman( 'run', matlabbatch );
      catch errMsg
        disp(errMsg); % probably contrast was not estimable due to it being based on not estimable betas
        figure; bar([conVec; est]'); legend({'Contrast','Estimable regressors'},'location','SouthWest')
        keyboard
      end
    end
    fclose('all'); % to avoid "Too many files open" error
  end
end


