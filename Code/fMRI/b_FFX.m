%% FFX batch code for Responsibility project. Can run "main" (conventional) and computational models.
%
% J Schultz 2019-2024

clear, %close all
warning('off','MATLAB:dispatcher:nameConflict')
spm12

%% ----------- Main settings ----------------------------------------------

% ----- Determine task(s)
% tasks = {'checkRealignmentParameters'};
% tasks = {'checkNtrials'};
% tasks = {'calculateModelFit'};
% tasks = {'FFXspecAndEst'};
% tasks = {'contrastCheck'};
tasks = {'contrasts'};
% tasks = {'FFXspecAndEst','checkNtrials'};
% tasks = {'FFXspecAndEst','contrasts'};

%% ---- Code and data directories
[~, computerID] = system('hostname');
if strfind(computerID,'Mac-Pro.local') % then it's my Mac Pro at CENs
  mainFolderData = '/Volumes/LaCie_RAID_A/Documents/exp/SODEC/';
else % Laptop; assumes it's connected to WD Elements external HD
  mainFolderData = '/Volumes/Elements/MasterMariaTom/SODEC/';
end
codeDir = '/Users/johannesschultz/sciebo/CURRENT WORK/Responsibility/Code/';
addpath(genpath(codeDir))

% ---- Mask for model fit calculations
anatMask = '/Users/johannesschultz/sciebo/NIFTIbrainMasks/Anatomy/StriatumBilat.nii';

%% ------------ Subjects --------------------------------------------------
% participants{1} = [ 44 55 62 63 65 69 70 71 74 77 78 79 81 83 84 85 87 90 92 93 94 95 96 98 101 102 103 104 106 110 120 121 124 125 126 129 136 155 166 171 173 174 175 176 177 ]; % all participants
% Subjects to do:
participants{1} = [ 44 55 62 63 65 69 70 71 74 77    79 81 83 84 85    90 92 93 94 95 96 98 101 102 103 104 106 110 120     124 125 126 129 136 155 166 171 173 174 175 ]; % all participants except those without data in OFC (87 and 176), people with only 1 session (121) or aborted Sess 1 (78), too few "risky" choices (87), missing risky outcome other pos social and self neg solo in Sess 1 (177). Still in are 62 (only 1 trial of outcome other neg partner in Sess 1), 65 and 120, missing outcome self neg solo in both sess (65) and Sess 1 (120)
participants{1} = [ 65 120 ]; % test participant without outcome risky neg n-soc trials
Ns = sum(cellfun(@length,participants)); % N subjects

%% --------- FFX model options --------------------------------------------
whichFFXmodel = 1; % 1 is main model, 2 is computational model.

switch whichFFXmodel % can do main FFX and comp models. gPPI done in another script.
  
  case 1    % The main model
    trialTypeCode = 'bb_onsetMaker';
    startTimeAndSetupValuesFile = [codeDir 'ResponsibilityExp_startTime.txt'];
    FFXfolderName = {'Responsibility/FFX/main/'};
    NtrialTypes = 19;
    paraMod = 0; NparamMod = 0;
    compRegressors = 0;
    % Regressors: 3 options, 6 decisions, 9 outcomes, 1 happiness, then 6 realignment params, total = 25 per run. We will only specify for the 1st session, code will complete for both.
    regressors = {...
      1, 'Options social';...
      2, 'Options partner';...
      3, 'Options solo';...
      4, 'Decision safe social';...
      5, 'Decision safe partner';...
      6, 'Decision safe solo';...
      7, 'Decision risky social';...
      8, 'Decision risky partner';...
      9, 'Decision risky solo';...
      10, 'Outcome safe social';...
      11, 'Outcome safe partner';...
      12, 'Outcome safe solo';...
      13, 'Outcome other pos social';...
      14, 'Outcome other pos partner';...
      15, 'Outcome self pos solo';...
      16, 'Outcome other neg social';...
      17, 'Outcome other neg partner';...
      18, 'Outcome self neg solo';...
      19, 'happiness';...
      };
    NRegressors = size(regressors,1);
    
  case 2    % Reduced model with onsets only for happiness ratings, and comp regressors
    trialTypeCode = 'ResponsibilityExpOnsetMaker';
    startTimeAndSetupValuesFile = [codeDir 'ResponsibilityExp_startTime.txt'];
    FFXfolderName = {'Responsibility/FFX/comp/'};
    NtrialTypes = 1;
    paraMod = 0; NparamMod = 0;
    compRegressors = 1;  regressorFolder = 'fMRIregressorsOtherRPEactPassModel';
    % Regressors: 1 happiness ratings & 5 comp regressors = 6, + 6 realignment parameters.
    regressors = {...
      1, 'happiness';...
      2, 'CR';...
      3, 'EV';...
      4, 'sRPE';...
      5, 'pRPEsocial';...
      6, 'pRPEpartner';...
      };
    NRegressors = size(regressors,1);
    
end % the 2 models

durationsAll = 0;
printCons = 0;   conNo = 1; % print contrast results to ps file?
HRFderivatives = [0 0];
TR = 2.5;
imgType = 's6wuf';

%% ------------ Run on Terminal and batch engagement options --------------
processDuration = NaN;
nojvm = 0; % set to 1 for -nojvm start
runBatch = 1; % Set to 0 to not execute SPM batch

%%  ---------- Contrasts --------------------------------------------------
deleteCons = 1; % deletes all already calculated contrasts

% T-test over all conditions:
psfileNewNameRoot = 'SPMresults ';
conType = 'T';

Ncond = NtrialTypes; Nbf = sum(HRFderivatives)+1; Nruns = 4; NregPerRun = Ncond*Nbf*(paraMod+2)+6;

switch whichFFXmodel % can do main FFX and comp models. gPPI done in another script.
  % Contrasts for RFX are specified here for 1 run; they get expanded for more runs below.
  
  case 1    % The main model
    
    % I only want the contrasts for decisions and outcomes, so numbers 4 to 18
    possibleSingleConditionContrasts = [eye(NRegressors) zeros(NRegressors,6)]; % all the theoretically possible single-condition contrasts
    cons = {}; i = 1;
    for c = 4:18
      cons{i,1} = regressors{c,2}; % names are filled in later
      cons{i,2} = possibleSingleConditionContrasts(c,:);
      i = i+1;
    end
    removeMissingConditionsFromContrast = 1; % removes regressors with missing onsets from the contrast specification
    
    % 4 specific contrasts: Decisions risky>safe and social>solo, outcome risky>safe, guilt effect
    %       cons{1,1} = 'Decision risky > safe';   con = zeros(1,19+6); con(7:9) = 1; con(4:6) = -1; cons{1,2} = con;
    %       cons{2,1} = 'Decision social > solo';  con = zeros(1,19+6); con([4 7]) = 1; con([6 9]) = -1; cons{2,2} = con;
    %       cons{3,1} = 'Outcome risky > safe';    con = zeros(1,19+6); con(13:18) = 1; con(10:12) = -2; cons{3,2} = con;
    % cons{i,1} = 'Guilt effect - Outcome other neg social > partner';  con = zeros(1,25); con(16) = 1; con(17) = -1; cons{i,2} = con;
    
  case 2    % Reduced model with onsets only for happiness ratings, and comp regressors
    
    % I only want the contrasts for the computational model regressors, so numbers 2 to 6 
    possibleSingleConditionContrasts = [eye(NRegressors) zeros(NRegressors,6)]; % all the theoretically possible single-condition contrasts
    cons = {}; i = 1;
    for c = 2:6
      cons{i,1} = regressors{c,2}; % names are filled in later
      cons{i,2} = possibleSingleConditionContrasts(c,:);
      i = i+1;
    end
    removeMissingConditionsFromContrast = 0;
    
    % 2 specific contrasts
    %     cons = {};
    %     cons{1,1} = 'CR+EV';   con = zeros(1,6+6); con(2:3) = 1; cons{1,2} = con;
    %     cons{2,1} = 'pRPEsocial - pRPEpartner';   con = zeros(1,6+6); con(5) = 1; con(6) = -1; cons{2,2} = con;
    
end % the 2 models

consToDo = 1:size(cons,1);

%% Settings unlikely to change
% order of the runs, for each subject, with name used for fMRI data and for behavioural file:
runNames = {...
  'rutledge1','run1';...
  'rutledge2','run2';...
  };
Nses = size(runNames,1);

%% Setup finished, get to work
for tt = 1:length(tasks)
  reportFile = [codeDir '.txt'];
  crashedList = [];
  doWhat = tasks{tt};
  disp('_______________________________________________')
  disp('|                                             |')
  disp(['| NOW DOING: ' doWhat ])
  disp('|_____________________________________________|')
  if runBatch
    spm( 'defaults', 'fmri' );
    spm_jobman( 'initcfg' );
    if ~nojvm;
      f_hdle = figure(99); f_hdle.Position = [130 877 560 420]; f_hdle.Color = [1 1 1]; g_hdle = gca; g_hdle.Visible='off'; drawnow; pause(.1)
    end
  end
  
  counter = 0;
  if ~nojvm; progressbar(0,0,doWhat); end
  for person_i = 1 : length( participants{1} )
    t1 = clock;
    % determine person name
    personNo = participants{1}(person_i);
    personName = ['SODEC_FMRI_' num2str( personNo ) ];
    fprintf('Doing participant: %s\n',personName)
    counter = counter + 1;
    
    try % so that if analysis for 1 subject crashes, the whole batch analysis does not stop
      % report
      if runBatch
        fid = fopen( reportFile, 'a' );
        string = sprintf('task %d/%d\ndoing %s\nsubj %d/%d\n%s\nLast model took %d seconds',...
          tt,length(tasks),doWhat,person_i,Ns,personName,processDuration);
        if ~nojvm;
          figure(99); cla; text(.01,.5,string,'fontsize',36,'color',[0 0 1],'interpreter','none'); g_hdle = gca; g_hdle.Visible='off'; drawnow
        end
      end
      
      % ===================== what to do? =======================
      switch doWhat
        case 'checkRealignmentParameters'
          figure;
          for ses = 1:Nses
            runNameFull = [mainFolderData '/fMRIdataComplete' filesep personName filesep runNames{ses,1}];
            multiregFile = dir( [runNameFull '/rp*.txt' ] );
            if length(multiregFile)>0
              sessMultiregFile = [runNameFull filesep multiregFile.name];
            else
              disp('No realignment parameters file found!')
            end
            if Nses > 1
              subplot(round(sqrt(Nses)),round(sqrt(Nses)),ses)
            end
            rp = dlmread(sessMultiregFile);
            plot(rp)
            axis tight
            set(gca,'ylim',[-4 4])
            title([runNames{ses,1}])
          end
          if Nses>1
            suptitle([personName])
          else
            title([personName])
          end
          if isempty(find(abs(rp)>3))
            printfig(['realignmentParameters_' personName])
          else
            printfig(['realignmentParameters_' personName '_over3'])
          end
          close
          
        case 'calculateModelFit'
          %           error('deactivated')
          % multiply residual model mask and anatomical mask images
          % with each other and average all voxels
          maskImg = [FFXdirName 'mask.nii'];
          ResMS = [FFXdirName 'ResMS.nii'];
          clear matlabbatch
          matlabbatch{1}.spm.util.imcalc.input = {maskImg; ResMS; anatMask};
          matlabbatch{1}.spm.util.imcalc.output = 'temp';
          matlabbatch{1}.spm.util.imcalc.outdir = {codeDir};
          matlabbatch{1}.spm.util.imcalc.expression = 'i1.*i2.*i3';
          matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
          matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
          matlabbatch{1}.spm.util.imcalc.options.mask = 0;
          matlabbatch{1}.spm.util.imcalc.options.interp = 1;
          matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
          spm_jobman( 'run', matlabbatch );
          img = spm_read_vols(spm_vol([codeDir 'temp.nii']));
          modelFit(ses,ISI) = mean(img(find(img)));
          
        case 'checkNtrials'
          % ------------ for each run type: ------------------
          for ses = 1:Nses
            clear scanInput sessMultiFile sessMultiregFile
            
            % Get behavioural data
            behavFile = dir([mainFolderData '/fMRIdataComplete/' personName filesep '*Decision_Rutledge*' runNames{ses,2} '*.mat']);
            fname = [behavFile.folder filesep behavFile.name];
            initialISIandSetupValues = dlmread(startTimeAndSetupValuesFile);
            idx = find(initialISIandSetupValues(:,1)==personNo);
            ISI = initialISIandSetupValues(idx,ses+1);
            if isempty(ISI); ISI = NaN; end
            ISIs(ses) = ISI;
            eval(['[onsets,durations,trial_type_names,param_mods,missingOnsets(person_i,ses,:)] = ' trialTypeCode '(fname,ISI);']);
            disp(['N trials per condition: ' sprintf('%d ',cellfun(@length,onsets))])
          end
          dlmwrite(['missingOnsets_' trialTypeCode '.txt'],double([participants{1}' squeeze(missingOnsets(:,1,:)+missingOnsets(:,2,:))]));
          
        case 'FFXspecAndEst'
          
          % make FFX folder
          FFXdirName = [mainFolderData filesep FFXfolderName{1} filesep personName filesep];
          if ~exist( FFXdirName, 'dir' ); mkdir(FFXdirName); fprintf(fid, 'Erstelle Ordner %s/r/n', FFXdirName ); end
          fprintf('FFX dir name: %s\n',FFXdirName)
          
          % ------------ for each run type: ------------------
          clear matlabbatch sessionData regs
          for ses = 1:Nses
            clear scanInput sessMultiFile sessMultiregFile
            
            % Get scans
            scanFolder = [mainFolderData '/fMRIdataComplete' filesep personName filesep runNames{ses,1} filesep];
            fprintf('fMRI data folder: %s\n',scanFolder)
            scanListing = dir( [scanFolder imgType '*.nii'] );
            scanInput = cell( length( scanListing ), 1 );
            for scanFile = 1 : length( scanListing )
              scanInput{ scanFile } = [scanFolder scanListing( scanFile ).name ',1'];
            end
            scans = sort( scanInput );
            sessionData{ses,1} = scans;
            
            % Get behavioural data for this experiment
            behavFile = dir([mainFolderData '/fMRIdataComplete' filesep personName filesep '*Decision_Rutledge*' runNames{ses,2} '*.mat']);
            fname = [behavFile.folder filesep behavFile.name];
            fprintf('Behaviour file: %s\n',fname)
            initialISIandSetupValues = dlmread(startTimeAndSetupValuesFile);
            idx = find(initialISIandSetupValues(:,1)==personNo);
            ISI = initialISIandSetupValues(idx,ses+1);
            if isempty(ISI); ISI = NaN; end
            ISIs(ses) = ISI;
            eval(['[onsets,durations,trial_type_names,param_mods,missingOnsets(person_i,ses,:)] = ' trialTypeCode '(fname,ISI);']);
            disp(['N trials per condition: ' sprintf('%d ',cellfun(@length,onsets))]) %#ok<*DSPS>
            
            % if no clean initialISI, set missingOnsets for all trialtypes to 1
            if isnan(ISI)
              missingOnsets(person_i,ses,:) = 1;
            end
            
            if whichFFXmodel == 2 % keep only the happiness ratings as conditions
              trial_type_names = trial_type_names(19);
              onsets = onsets(19);
              durations = durations(19);
              param_mods = param_mods(19);
            end
            
            sessionData{ses,2} = trial_type_names;
            sessionData{ses,3} = onsets;
            sessionData{ses,4} = durations;
            
            % Get realignment parameters
            runNameFull = [mainFolderData '/fMRIdataComplete' filesep personName filesep runNames{ses,1}];
            multiregFile = dir( [runNameFull '/rp*.txt' ] );
            if length(multiregFile)>0
              sessMultiregFile = [runNameFull filesep multiregFile.name];
            end
            sessionData{ses,5} = sessMultiregFile;
            
            % Get parametric modulation variables
            sessionData{ses,6} = param_mods;
            
            % Get computational regressors
            if compRegressors
              % load parametric modulator regressors
              regs = dlmread([codeDir regressorFolder '/fMRIregressors_' num2str(personNo) '_' num2str(ses) '.txt']);
              sessionData{ses,7} = regs;
            end
            
          end % session
          
          if length(sessionData(1,:)) ~= 6; warning('Input data missing!');
            sessionData
          end
          
          if ~runBatch
            sessionData
          else
            % clean the FFX directory
            owd = pwd; cd(FFXdirName); delete *; cd(owd);
          end
          
          % throw away session information for session with no ISI
          sessionData = sessionData(find(~isnan(ISIs)),:);
          
          % specify multi-session model
          matlabbatch{1}.spm.stats.fmri_spec.dir = { FFXdirName };
          matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
          matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
          matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
          matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
          for ses = 1:size(sessionData,1) % add session-specific infos
            matlabbatch{1}.spm.stats.fmri_spec.sess(ses).scans = sessionData{ses,1};
            % load(sessionData{ses,3}) % loads onsets, durations and difficulties (-> parametric modulators) into workspace
            % matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond = struct('name', names, 'onset', onsets, 'duration', durations, 'tmod', {}, 'pmod', struct('name', {'difficulty'}, 'param', difficulties, 'poly', {1}));
            for c = 1:length(sessionData{ses,2}) % for each condition
              matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).name = sessionData{ses,2}{c};
              matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).onset = sessionData{ses,3}{c};
              %               matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).duration = durationsAll(tt) * ones(length(sessionData{ses,3}{c}),1);% durations{c};
              matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).duration = sessionData{ses,4}{c};
              matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).tmod = 0;
              if paraMod && ~isempty(sessionData{ses,6}{c}) && condWithParamod(c)
                matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).pmod = struct('name', {'Val'}, 'param', sessionData{ses,6}{c}, 'poly', {1});
              else
                matlabbatch{1}.spm.stats.fmri_spec.sess(ses).cond(c).pmod = struct('name', {}, 'Val', {});
              end
            end
            if compRegressors
              for rrr = 1:size(sessionData{ses,7},2)
                NtimePoints = length(sessionData{ses,1}); % this is Nscans
                dif = size(sessionData{ses,7},1) - NtimePoints;
                if dif < 0; sessionData{ses,7} = [sessionData{ses,7}; zeros(ceil(-dif),5)]; end
                matlabbatch{1}.spm.stats.fmri_spec.sess(ses).regress(rrr).name = regressors{rrr+1,1};
                matlabbatch{1}.spm.stats.fmri_spec.sess(ses).regress(rrr).val = sessionData{ses,7}(1:NtimePoints,rrr);
              end
            else
              matlabbatch{1}.spm.stats.fmri_spec.sess(ses).regress = struct('name', {}, 'val', {});
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess(ses).multi_reg = sessionData(ses,5);
            matlabbatch{1}.spm.stats.fmri_spec.sess(ses).hpf = 128;
          end
          matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
          matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = HRFderivatives;
          matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
          matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
          matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
          matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
          
          % estimate model
          matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = { [FFXdirName 'SPM.mat'] };
          matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
          save([FFXdirName '/sessionData'],'sessionData')
          save([FFXdirName '/matlabbatch'],'matlabbatch')
          if runBatch
            spm_jobman( 'run', matlabbatch );
          end
          
          % ================== CONTRASTS ===========================
        case 'contrastCheck'
          % check whether the contrasts are estimable
          FFXdirName = [mainFolderData filesep FFXfolderName{1} filesep personName filesep];
          
          clear matlabbatch;
          matlabbatch{1}.spm.stats.con.spmmat = { [FFXdirName 'SPM.mat'] };
          matlabbatch{1}.spm.stats.con.delete = deleteCons;
          
          % first, determine for each contrast if the necessary conditions
          % have onsets
          initialISIandSetupValues = dlmread(startTimeAndSetupValuesFile);
          idx = find(initialISIandSetupValues(:,1)==personNo);
          ISIs = initialISIandSetupValues(idx,2:end);
          clear missingOnsetsForCon
          for ses = 1:Nses
            behavFile = dir([mainFolderData '/fMRIdataComplete' filesep personName filesep '*Decision_Rutledge*' runNames{ses,2} '*.mat']);
            fname = [behavFile.folder filesep behavFile.name];
            eval(['[~,~,trial_type_names,~,missingOnsetsForCon(ses,1:NtrialTypes)] = ' trialTypeCode '(fname,NaN);']);
            
            % if no clean initialISI, set missingOnsets for all trialtypes to 1
            if isnan(ISIs(ses))
              missingOnsetsForCon(ses,:) = 1;
            end
          end
          
          % prepare contrast names
          try load([FFXdirName 'SPM.mat'])
            conNames = char(SPM.xX.name(1:size(cons,1)));
            for c = 1:size(cons,1)
              con1run = cons{c,2}; % get contrast values
              if Nses > 1
                cons{c,1} = deblank(conNames(c,7:end)); % get contrast names
              else
                cons{c,1} = deblank(conNames(c,:)); % get contrast names
              end
              
              % remove sessions without these trials
              runMultiplicator = ones(1,Nses);
              if removeMissingConditionsFromContrast
                for ses = 1:Nses
                  if sum(missingOnsetsForCon(ses,find(con1run(1:NtrialTypes+(NparamMod+1)))))
                    runMultiplicator(ses) = 0;
                  end
                end
                runMultiplicator = runMultiplicator(find(~isnan(ISIs)));
              end
              
              % build contrast vector
              conVec = [kron(runMultiplicator,con1run) zeros(1,length(runMultiplicator))];
              
              % remove missing conditions
              if removeMissingConditionsFromContrast
                mm = [makerow([missingOnsetsForCon zeros(Nses,6)]') zeros(1,Nses)];
                conVec = conVec(~mm);
              end
              
              % check if contrast is estimable
              [~,~,emsg{person_i,c},imsg] = spm_conman('ParseCon',conVec,SPM.xX.xKXs,conType);
              %             if ~isempty(emsg{person_i,c})
              %               disp(emsg{person_i,c});
              %             end
            end % contrast
            %           keyboard
          catch
            disp([FFXdirName 'SPM.mat cannot be loaded'])
          end
        case 'contrasts'
          % set the contrasts for each subject depending on the conditions available
          FFXdirName = [mainFolderData filesep FFXfolderName{1} filesep personName filesep];
          
          clear matlabbatch;
          matlabbatch{1}.spm.stats.con.spmmat = { [FFXdirName 'SPM.mat'] };
          matlabbatch{1}.spm.stats.con.delete = deleteCons;
          
          % first, determine for each contrast if the necessary conditions have onsets
          initialISIandSetupValues = dlmread(startTimeAndSetupValuesFile);
          idx = find(initialISIandSetupValues(:,1)==personNo);
          ISIs = initialISIandSetupValues(idx,2:end);
          clear missingOnsetsForCon
          for ses = 1:Nses
            behavFile = dir([mainFolderData '/fMRIdataComplete' filesep personName filesep '*Decision_Rutledge*' runNames{ses,2} '*.mat']);
            fname = [behavFile.folder filesep behavFile.name];
            if whichFFXmodel == 1 % in comp model, the happiness rating is the only condition
              eval(['[~,~,trial_type_names,~,missingOnsetsForCon(ses,1:NtrialTypes)] = ' trialTypeCode '(fname,ISIs(ses));']);
            end
            % if no clean initialISI, set missingOnsets for all trialtypes to 1
            if isnan(ISIs(ses))
              missingOnsetsForCon(ses,:) = 1;
            end
          end
          
          for c = 1:size(cons,1)
            con1run = cons{c,2}; % get contrast values
            conName = deblank(cons{c,1}); % get contrast names
            
            % remove sessions without these trials
            runMultiplicator = ones(1,Nses);
            if removeMissingConditionsFromContrast
              for ses = 1:Nses
                if sum(missingOnsetsForCon(ses,find(con1run(1:NtrialTypes+NparamMod))))
                  runMultiplicator(ses) = 0;
                end
              end
              runMultiplicator = runMultiplicator(find(~isnan(ISIs)));
            end
            
            % build contrast vector
            conVec = [kron(runMultiplicator,con1run) zeros(1,length(runMultiplicator))];
            
            % remove missing conditions from the contrast
            if removeMissingConditionsFromContrast
              keep = 1-[makerow([missingOnsetsForCon zeros(Nses,6)]') zeros(1,Nses)]; % one if regressor has no missing onsets
              conVec = conVec.*keep; % regressors with missing trials are weighted with 0 in any contrast, thus deselecting them if they would have been included
            end
            
            if sum(abs(conVec)) % at least one non-zero value
              disp(['Running contrast: ' cons{c,1}])
              if c > 1;  matlabbatch{1}.spm.stats.con.delete = 0;  end
              switch conType
                case 'F'
                  matlabbatch{1}.spm.stats.con.consess{1}.fcon.name = conName;
                  matlabbatch{1}.spm.stats.con.consess{1}.fcon.convec = conVec;
                  matlabbatch{1}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
                case 'T'
                  matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = cons{c,1};
                  matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = conVec;
                  matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
              end
              if printCons
                matlabbatch{2}.spm.stats.results.spmmat = { [FFXdirName 'SPM.mat'] };
                matlabbatch{2}.spm.stats.results.conspec.titlestr = '';
                matlabbatch{2}.spm.stats.results.conspec.contrasts = conNo(c);
                matlabbatch{2}.spm.stats.results.conspec.threshdesc = 'none';
                matlabbatch{2}.spm.stats.results.conspec.thresh = 0.001;
                matlabbatch{2}.spm.stats.results.conspec.extent = 5;
                matlabbatch{2}.spm.stats.results.conspec.mask = struct('contrasts', {}, 'thresh', {}, 'mtype', {});
                matlabbatch{2}.spm.stats.results.units = 1;
                matlabbatch{2}.spm.stats.results.print = true;
              end
              if runBatch
                spm_jobman( 'run', matlabbatch );
              end
              if printCons
                psfile = dir([FFXdirName filesep '*.' outputFileType]);
                seps = findstr(FFXdirName,filesep);
                newDir = FFXdirName;
                newDir(seps([5 6])) = '_';
                copyfile([FFXdirName filesep psfile(1).name],[newDir psfileNewNameRoot cons{c,1} '.' outputFileType])
                delete(psfile(1).name)
              end
            else
              disp(['Contrast not estimated because of missing trials: ' cons{c,1}])
            end % not empty contrast
          end % each con
      end % doWhat
      % post batch cleanup
      switch doWhat
        case 'contrasts'
          % rename ps output file
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      t2 = clock;
      processDuration = floor(etime( t2,t1));
      if runBatch
        fprintf(fid, [datetime_js 'Duration for participant ' personName ' in seconds: %d\r\n'], processDuration );
        fprintf(fid, '\r\n' );
      end
      report{personNo} = [personName ' ran fine'];
      
    catch errorMsg
      crashedList = [crashedList personNo];
      reportedErrorMessages{personNo} = errorMsg;
      report{personNo} = [personName ' crashed: ' errorMsg.message];
      disp(errorMsg.stack(end))
      %                 keyboard
    end % try/catch
    % post batch cleanup
    switch doWhat
      case 'calculateModelFit'
        figure('name',['ModelFit ' personName ', ' runNames{ses}]);
        plot(initialISI,modelFit(ses,initialISI),'o-')
        xlabel('ISIs [s]')
        ylabel('mean residuals in Area V1')
        [~,idx] = min(modelFit(ses,initialISI));
        bestFittingISI(personNo,ses) = initialISI(idx);
        owd = pwd; cd(codeDir)
        printfig
        close
        cd(owd)
    end
    if ~nojvm; progressbar(person_i / length( participants{1} ),0,[doWhat ': ' personName]); end
    
  end % loop over participants
  
  switch doWhat
    case 'checkNtrials'
      figure('position',[1000         630         876         708]);
      imagesc(squeeze(missingOnsets(:,1,:)+missingOnsets(:,2,:))); colorbar
      set(gca,'XTick',1:length(trial_type_names)); xlabel_oblique(trial_type_names,40,12,1)
      set(gca,'ytick',1:size(missingOnsets,1),'yticklabel',participants{1})
      set(gca,'position',[0.1184    0.2762    0.7058    0.6488])
      suptitle('N sessions without onsets for each participant')
    case 'contrastCheck'
      %       disp(emsg)
      badOnes = ~cellfun(@isempty,emsg);
      figure; imagesc(badOnes)
      set(gca,'XTick',1:size(cons,1)); xlabel_oblique(cons(:,1),40,12,1)
      set(gca,'position',[0.1300    0.3857    0.7750    0.5393])
  end
  
  
  disp('....----====|||||====----....')
  if ~isempty(crashedList)
    disp(['Crashed: ' num2str(crashedList)])
  else
    disp('All jobs ran without crash')
  end
end % tasks

if runBatch
  fprintf(fid, 'Done!' );
  fclose( fid );
end
