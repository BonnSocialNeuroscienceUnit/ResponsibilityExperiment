%% SPM12 RFX batch code for Responsibility project. Can run "main" (conventional) and computational models.
%
% J Schultz 2019-2024

clear
baseDir = '/Volumes/LaCie_RAID_A/Documents/exp/SODEC/Responsibility/';
FFXbaseDir = [baseDir 'FFX/'];
RFXbaseDir = [baseDir 'RFX/FullFact/']; % Decision and Outcome are in subfolders of this

% command flags
estimate = 1; % if 0, does not estimate the model
isPPI = 0; % set to 1 to run connectivity analyses
whichRFX = 'main/decision';   % main FFX model, decision regressors
% whichRFX = 'main/outcome';    % main FFX model, outcome regressors
% whichRFX = 'comp';            % model with computational regressors
whichJob = 'run'; %'interactive'; % if "interactive", runs through the spm batch system, allows to inspect the batch set up

%% Subject info
subjPrefix = 'SODEC_FMRI_';
% for main FFX decisions and comp model, I can take all the following participants:
subjNumbers = [ 301 302 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 347 334 335 336 337 338 339 340 ];
% for main FFX outcome, I have to remove participant 306 who had no outcome self neg solo in both sessions:
subjNumbers = [ 301 302 304 305     307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 347 334 335 336 337 338 339 340 ];

%% PPI information
if isPPI
  PPInames = {...
    'PPI_insLeft_outcRiskySafe',...
    'PPI_STSleft_sRPEsocPartner',...
    };
else
  PPInames = {''};
end

%% Input: contrast names in the FFX model that this RFX is based on, and factors for the RFX
switch whichRFX
  case 'main/decision' % also for PPIs
    % Which model?
    FFXname = 'main'; % also for PPIs, PPI name gets added

    % Contrast names and numbers
    FFXcontrastNumbers = [1 3 4 6]; % no need to examine partner decisions, as the participant didn't make them!
    FFXcontrastNames = {...
      'Decision safe social'...
      'Decision safe solo'...
      'Decision risky social'...
      'Decision risky solo'...
      };

    % Factors
    factNames = {'SafeRisky','SocialSolo'};
    levels = [...
      1 1;...
      1 2;...
      2 1;...
      2 2;...
      ];

    % Explicit mask
    explicitMask = {''};
    
    % Contrasts to estimate in the RFX (if I don't delete the automatically esimated ones, these are added at the end)
    RFXcontrasts(1).name  = 'Safe>Risky';
    RFXcontrasts(1).c     = [1 1 -1 -1]';
    RFXcontrasts(1).STAT  = 'T';
    RFXcontrasts(2).name  = 'Risky>Safe';
    RFXcontrasts(2).c     = [-1 -1 1 1]';
    RFXcontrasts(2).STAT  = 'T';
    RFXcontrasts(3).name  = 'Social>Solo';
    RFXcontrasts(3).c     = [1 -1 1 -1]';
    RFXcontrasts(3).STAT  = 'T';
    RFXcontrasts(4).name  = 'Solo>Social';
    RFXcontrasts(4).c     = [-1 1 -1 1]';
    RFXcontrasts(4).STAT  = 'T';
    RFXcontrasts(5).name  = 'Risky>Safe X Social>Solo';
    RFXcontrasts(5).c     = [-1 1 1 -1]'; % the result of [-1 -1 1 1].*[1 -1 1 -1]
    RFXcontrasts(5).STAT  = 'T';
    RFXcontrasts(6).name  = 'Risky>Safe X Solo>Social';
    RFXcontrasts(6).c     = [1 -1 -1 1]'; % the result of [-1 -1 1 1].*[-1 1 -1 1]
    RFXcontrasts(6).STAT  = 'T';
    
    % For display: contrast and thresholding information
    selectedContrast = 1;
    extentThresh = 100; % 100 to get only the FWE-surviving results
    SVCimage = '';
    
  case 'main/outcome'
    % Which model?
    FFXname = 'main';

    % Contrast names and numbers
    FFXcontrastNumbers = 7:15;
    FFXcontrastNames = {...
      'Outcome safe social';...
      'Outcome safe partner';...
      'Outcome safe solo';...
      'Outcome other pos social';...
      'Outcome other pos partner';...
      'Outcome self pos solo';...
      'Outcome other neg social';...
      'Outcome other neg partner';...
      'Outcome self neg solo';...
      };

    % Factors
    factNames = {'Safe_RiskyPos_RiskyNeg','SocialPartnerSolo'};
    levels = [...
      1 1;...
      1 2;...
      1 3;...
      2 1;...
      2 2;...
      2 3;...
      3 1;...
      3 2;...
      3 3;...
      ];
    
    % Explicit mask
    explicitMask = {''};
    
    % Contrasts to estimate in the RFX (if I don't delete the automatically esimated ones, these are added at the end)
    RFXcontrasts(1).name  = 'Risky>Safe';
    RFXcontrasts(1).c     = [-2 -2 -2 1 1 1 1 1 1]';
    RFXcontrasts(1).STAT  = 'T';
    RFXcontrasts(2).name  = 'NegSocial>NegPartner'; % gives weak results here
    RFXcontrasts(2).c     = [0 0 0 0 0 0 1 -1 0]';
    RFXcontrasts(2).STAT  = 'T';
    
    % For display: contrast and thresholding information
    selectedContrast = 1;
    extentThresh = 130; % to get only the FWE-surviving results
    SVCimage = '';

  case 'comp'
    % Which model?
    FFXname = 'comp';

    % Contrast names and numbers
    FFXcontrastNumbers = 1:5;
    FFXcontrastNames = {...
      'CR';...
      'EV';...
      'sRPE';...
      'pRPEsocial';...
      'pRPEpartner';...
      };

    % Factors
    factNames = {'CompRegs'};
    levels = [...
      1;...
      2;...
      3;...
      4;...
      5;...
      ];

    % Explicit mask
    explicitMask = {''};

    % Contrasts to estimate in the RFX (if I don't delete the automatically esimated ones, these are added at the end)
    RFXcontrasts(1).name  = 'CR+EV';
    RFXcontrasts(1).c     = [1 1 0 0 0]';
    RFXcontrasts(1).STAT  = 'T';
    RFXcontrasts(2).name  = 'pRPEsocial>pRPEpartner';
    RFXcontrasts(2).c     = [0 0 0 1 -1]';
    RFXcontrasts(2).STAT  = 'T';

    % For display: contrast and thresholding information
    selectedContrast = 2;
    extentThresh = 70; % to get only the FWE-surviving results
    SVCimage = [baseDir 'RFX/FullFact/main/outcome/risky>safe_0p001u_k130.nii']; % for contrast 2, but gives different results than using the GUI. No idea why.
    SVCimage = ''; % Use this for all contrasts

end

%% setup batch
if isempty(which('spm')) || ~strcmpi(spm('ver'),'SPM12'); spm12; spm( 'defaults', 'fmri' ); end
warning('OFF', 'MATLAB:dispatcher:nameConflict')
for pp = 1:length(PPInames) % one RFX per PPI that we run
  
  %% setup RFX model
  spm_jobman( 'initcfg' );
  spm('defaults','FMRI')
  clear matlabbatch
  if isPPI
    RFXdir = [fullfile(RFXbaseDir,[whichRFX '_' PPInames{pp}]) filesep];
  else
    RFXdir = [fullfile(RFXbaseDir,whichRFX) filesep];
  end
  disp(['RFX directory name: ' RFXdir])
  if ~exist(RFXdir,'dir'); mkdir(RFXdir); end
  matlabbatch{1}.spm.stats.factorial_design.dir = {RFXdir};
  
  % gather input information and set up model
  factNLevels = max(levels);
  
  for f = 1:length(factNames)
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).name = factNames{f};
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).levels = factNLevels(f);
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).dept = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).variance = 1;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).gmsca = 0;
    matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(f).ancova = 0;
  end
  
  % Set up input data
  con = FFXcontrastNumbers;
  for j=1:numel(con) % = N cells
    thiscell = {};
    for i=1:numel(subjNumbers)
      if isPPI
        thiscell{i,1} = fullfile(FFXbaseDir,FFXname,sprintf('SODEC_FMRI_%d',subjNumbers(i)),PPInames{pp},sprintf('con_%04d.nii',con(j)));
      else
        thiscell{i,1} = fullfile(FFXbaseDir,FFXname,sprintf('SODEC_FMRI_%d',subjNumbers(i)),sprintf('con_%04d.nii',con(j)));
      end
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(j).levels = levels(j,:)';
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(j).scans = thiscell;
  end
  matlabbatch{1}.spm.stats.factorial_design.des.fd.contrasts = 1;
  matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
  matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
  matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
  matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
  matlabbatch{1}.spm.stats.factorial_design.masking.em = explicitMask;
  matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
  matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
  matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
  
  % ----- Estimate if desired -----------------------------------------------
  if estimate
    delete([RFXdir 'SPM.mat']) % so that spm_run_factorial_design does not ask through GUI for permission to overwrite previous analysis
    matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
    matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
  end
  
  % ----- Add contrasts of interest -----------------------------------------
  matlabbatch{3}.spm.stats.con.spmmat = {[RFXdir 'SPM.mat']};
  matlabbatch{3}.spm.stats.con.delete = 1; % even if I don't delete the automatically created contrasts, additional contrasts are added afterwards
  
  for cc = 1:length(RFXcontrasts)
    if RFXcontrasts(cc).STAT == 'T'
      matlabbatch{3}.spm.stats.con.consess{cc}.tcon.name = RFXcontrasts(cc).name;
      matlabbatch{3}.spm.stats.con.consess{cc}.tcon.weights = RFXcontrasts(cc).c;
      matlabbatch{3}.spm.stats.con.consess{cc}.tcon.sessrep = 'none';
    else % F test
      matlabbatch{3}.spm.stats.con.consess{cc}.fcon.name = RFXcontrasts(cc).name;
      matlabbatch{3}.spm.stats.con.consess{cc}.fcon.weights = RFXcontrasts(cc).c;
      matlabbatch{3}.spm.stats.con.consess{cc}.fcon.sessrep = 'none';
    end
  end

  % Pick a contrast and threshold
  matlabbatch{4}.spm.stats.results.spmmat(1) = cfg_dep('Contrast Manager: SPM.mat File', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
  matlabbatch{4}.spm.stats.results.conspec.titlestr = '';
  matlabbatch{4}.spm.stats.results.conspec.contrasts = selectedContrast;
  matlabbatch{4}.spm.stats.results.conspec.threshdesc = 'none';
  matlabbatch{4}.spm.stats.results.conspec.thresh = 0.001;
  matlabbatch{4}.spm.stats.results.conspec.extent = extentThresh;
  matlabbatch{4}.spm.stats.results.conspec.conjunction = 1;
  if isempty(SVCimage)
    matlabbatch{4}.spm.stats.results.conspec.mask.none = 1;
  else
    matlabbatch{4}.spm.stats.results.conspec.mask.image.name = {SVCimage};
    matlabbatch{4}.spm.stats.results.conspec.mask.image.mtype = 0;
  end
  matlabbatch{4}.spm.stats.results.units = 1;
  matlabbatch{4}.spm.stats.results.export{1}.ps = false;
  
  % ----- Execute -----------------------------------------------------------
  spm_jobman(whichJob,matlabbatch);
end
