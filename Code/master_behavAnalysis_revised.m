%% Script to load and fit behavioural data from the Responsibility project experiments

%% General setup
clear; close all; clc
owd = pwd;
colors = niceColorsDark(3); % the colors for the 3 experimental conditions
fontSize = 20;
set(0,'defaultAxesFontSize',fontSize)
addpath('fittingUtilityFunction')

%% get datafiles
% By hand:
% During fMRI data:
mfilepath = [fileparts(mfilename('fullpath')) filesep];
addpath(genpath(mfilepath))
% baseDir = '/Users/johannesschultz/sciebo/CURRENT WORK/Responsibility/';
baseDir = [fileparts(fileparts(mfilepath)) filesep]; % should go 1 level higher than Code

%% select which study(ies) to show
study = {'Behav','fMRI','both (-> paper)'};
% try selection = bttnChoiceDialog(study, 'Which dataset?', 3,'',[1 2 3]);
%   study = study{selection};
% end
study = 'Behav';
study = 'fMRI';
study = 'both (-> paper)';

switch study
  case 'Behav'
    paperFigure = 0;
    datasets{1} = 'Behav';
    pathnames{1} = [baseDir 'BehaviouralData/BehaviourOnly']; % Behav only data
  case 'fMRI'
    paperFigure = 0;
    datasets{1} = 'fMRI';
    pathnames{1} = [baseDir 'BehaviouralData/fMRIstudy']; % fMRI expt data
  case 'both (-> paper)'
    paperFigure = 1;
    datasets{1} = 'Behav';
    datasets{2} = 'fMRI';
    pathnames{1} = [baseDir 'BehaviouralData/BehaviourOnly']; % Behav only data
    pathnames{2} = [baseDir 'BehaviouralData/fMRIstudy']; % fMRI expt data
end

%% Definition of condition names
% conditionNames = {'You decide for both','Partner decides for both','You decide for yourself'};
conditionNames = {'Social','Partner','Solo'};
hap_typeNames_soc = {'high | high','high | low','low | high','low | low'}; % first participant, then partner
hap_typeNames_nonsoc = {'high','low'};

%% display options
happyFigClassicOrDotPlotFlag = 2; % 0 -> no fig; 1 -> meansembar; 2 -> dotPlot. Choose 1 if want sequentialFiguresForTalk to be 1
sequentialFiguresForTalk = 0;
useDiffHappinessForSimplePlot = 1; % uses changes in subsequent happiness ratings for the simple behavioural happiness plot

if paperFigure
  studyOffsets            = [0 .45 0 0;0 0 0 0];
  paperFig1 = figure('name','Figure 2 - Choices','position',[5 5 1600 840]);
  decPanelPos             = [.06 .1 .40 .3];
  riskPremPanelPos        = [.55 .1 .16 .3];
  rhoPanelPos             = [.79 .1 .16 .3];
  paperFig2 = figure('name','Figure 3 - Happiness','position',[5 5 1600 840]);
  happyRewardPanelPos1    = [.06 .13 .14 .26];
  happyRewardPanelPos2    = [.28 .13 .14 .26];
  happyPredPanelPos       = [.50 .13 .14 .26];
  happyLotteryPanelPos    = [.72 .13 .22 .26];
end

%% Happiness data analysis settings
% Zscore the happiness results separately for each run
zscoreSepEachRun = 'Yes';

% % For the non-computational analysis, use only the trials directly before the happiness ratings?
useAllTrialsForNonComp = 'no';

% Computational modeling of happiness: do it or not?
compModelling = 1;
showIndivmodelNameFitFigs = 0;
useTauForFMRIreg = 1;

% Settings for computational modeling a la Rutledge
% Use which modelName?
if compModelling
  addpath('TDmodels')
  compModelName = {'None','Compare','ParamRecovery','Basic','InequalSimple','GuiltEnvy','otherRPEactPass','Responsibility','ResponsibilityRedux'};
  try selection = bttnChoiceDialog(compModelName,'',1,'Use which model name for computational analysis?',[3 3]);
    compModelName = compModelName{selection};
  end
  if strcmpi(compModelName,'None'), compModelling = 0; end
else; compModelName = 'None';
end

whichTrials = 'all';

% Create fMRI regressors?
compToDoList = []; % nothing per default
createfMRIreg = 'No'; % not per default
if compModelling && ~strcmpi(compModelName,'Compare') && ~strcmpi(compModelName,'ParamRecovery') && ~strcmpi(compModelName,'None')
  createfMRIreg = {'Yes','No'};
  try selection = bttnChoiceDialog(createfMRIreg,'',1,'Create fMRI regressors?',[1 2]);
    createfMRIreg = createfMRIreg{selection};
  end
end

% Determine comp modelling to dos
if compModelling, compToDoList = 1; end % default for comp modelling
if compModelling && strcmpi(compModelName,'ParamRecovery'), compToDoList = 0; end % no real fitting needed, but loading of fitted data and simulations
if strcmpi(createfMRIreg,'Yes')
  compToDoList = [1 2]; % modelling AND make fMRI regressors
  spm12 % put spm functions in path
  initialISIandSetupValues = dlmread([mfilepath 'bin/ResponsibilityExp_startTime.txt']);
end

%% Decision data analysis settings

% Task: choice between safe option and lottery with 50% probability of high
% vs. low outcome, in 2 contexts: decide for self only, or for self and
% other.
% Aim of analysis: calculate the certainty equivalents (CEs) of the
% lotteries in the 2 choice conditions. CEs are defined as the amount of
% reward for which the participant was choice-indifferent with regards to
% said lotteries; the CE therefore indicates the subjective value of the
% lottery in the currency of the safe reward. Lottery CEs larger than said
% lotteries' objective expected value (EV) reflect risk-seeking behavior;
% CEs smaller than the lotteries' EV indicate risk aversion.
%
% Then, estimate the contribution of utility to the subjective
% values in both choice conditions.

% analyseDecisions = {'Yes','No'};
% try selection = bttnChoiceDialog(analyseDecisions,'',1,'Analyse decisions?',[1 2]);
%   analyseDecisions = analyseDecisions{selection};
% end
analyseDecisions               = 'Yes'; % general problem: I have few trials and people are not always consistent. See table decTab for a subject's decisions as function of values of risky and safe options
flag.fitBinnedEVdiffs          = 1;  % the first type of decision analysis I did, abandoned for better alternatives
flag.fitBinnedSVdiffs          = 1;  % same, but based on subjective value. Only works if I have CARA rhos for each trialtype, so need flag.findCARArhoInSearchSpace=1
flag.fitProbitFunction         = 1;  % fits the choices fine as function of EVrisky-Vsafe, but (1) is this a valid independent value? 2) I can't directly interpret the estimated parameters
flag.findCARArhoInSearchSpace  = 1;  % this is the best so far: finds the risk aversion parameter that best explains the choices.
rhoToTest                      = -.1:0.001:.1; % the values of rho to test
flag.checkRiskAversionConsist  = 1;  % Only works if fitProbitFunctionFlag and findCARArhoInSearchSpaceFlag are both = 1

binEdges = -50:20:50;

choiceLMMs = {}; % we will collect choice probit linear mixed models in here

trialTypeColors = niceColors70s(3);

% Exponential utility function for the CARA model
% NOTE: rho is a constant that represents the degree of risk preference
% (rho>0 for risk aversion, rho=0 for risk-neutrality, and rho<0 for
% risk-seeking). c is a commodity, i.e. something the decision-maker likes
u_CARA   = @(c, rho)   (rho ~= 0) .* (1 - exp(-rho .* c)) ./ (   rho + (rho == 0) .* 0.000001)  + ...
  (rho == 0) .*                  c;
% ---- Expected utility function, which is the utility times the probability
% c1 and c2 are the high and low amounts of the risky option
EU_CARA   = @(c1, c2, p1, rho)  p1 .* u_CARA(c1, rho) + (1 - p1) .* u_CARA(c2, rho);


%% ------------------------------------------------------------------------
%
%                       END SETTINGS, NOW GET DATA
%
% -------------------------------------------------------------------------

% which dataset to use has been defined above

for currentStudy = 1:length(datasets)
  pathname = pathnames{currentStudy};
  dataset = datasets{currentStudy};
  studyName = sprintf('Study %d', currentStudy);
  
  if strcmpi(study,'both (-> paper)') & paperFigure
    studyOffset = studyOffsets(currentStudy,:);
  end
  fprintf(1,' =========== Dataset %s ============ \n',dataset)
  clear R RR subjNumbTemp initialISIandSetup filesCell
  
  files = dir([pathname '/*.mat']); filesCell = {};
  for f = 1:length(files); filesCell{f} = files(f).name; end; files = filesCell;
  % if ~iscell(files); ff{1} = files; files = ff; clear ff; end % needed for Windows I think
  
  % sort or not?
  if length(files) > 1
    %   sortOrNot = {'sort','don''t sort'};
    %   try selection = bttnChoiceDialog(sortOrNot, 'Sort according to subject number or not?', '','',[1 2]);
    %     sortOrNot = sortOrNot{selection};
    %   end
    sortOrNot = 'sort';
    switch sortOrNot
      case 'sort'
        files = sortrows(char(files));
      case 'don''t sort'
        files = char(files);
    end
  else
    files = char(files);
  end
  % files = files(1:10,:);
  
  %% load data
  for f = 1:size(files,1)
    Raw = load([pathname filesep deblank(files(f,:))]);
    Raw.results(1).subjName = deblank(files(f,:));
    Raw.results(1).Ntrials = length(Raw.results);
    ii = findstr(Raw.results(1).subjName,'_sub');
    temp = Raw.results(1).subjName(ii+5:ii+7);
    if strcmpi(temp(1),'M') || strcmpi(temp(1),'W') % then it's from Julia's experiment
      temp2 = str2num(temp(2:end));
      if strcmpi(temp(1),'M')
        temp2 = 100 + temp2;
      elseif strcmpi(temp(1),'W')
        temp2 = 200 + temp2;
      end
      temp = temp2;
    end % then it's not Julia's experiment, we should have the number
    if isempty(temp)
      temp = Raw.results(1).subjName(ii(end-2)+1:ii(end-1)-1);
    end
    if ischar(temp); temp=str2num(temp); end
    for str_loop = 1:length(Raw.results)
      gl=[Raw.results.amountThisTrialSubject];
      Raw.results(str_loop).Agl_Abs=sum(gl(1:str_loop-1));
      riskchoices=Raw.results(str_loop).riskyOption;
      riskEV=(riskchoices*0.5);
      choicesEV=[riskEV Raw.results(str_loop).safeOption];
      minEV=min(choicesEV);
      Raw.results(str_loop).MinEV=minEV;
      if Raw.results(str_loop).play==1
        choiceEV=sum(riskchoices*0.5);
        Raw.results(str_loop).ChoiceEV=choiceEV;
        Raw.results(str_loop).Ch_diff=choiceEV-minEV;
        
      elseif Raw.results(str_loop).play==0
        choiceEV=Raw.results(str_loop).safeOption;
        Raw.results(str_loop).ChoiceEV=choiceEV;
        Raw.results(str_loop).Ch_diff=choiceEV-minEV;
      end
      Raw.results(str_loop).PE=(gl(str_loop)-choiceEV);
    end
    subjNumbers(f) = temp;
    R{f} = Raw.results;
    if strcmpi(zscoreSepEachRun,'Yes')
      % zscore happiness (separately for each run)
      happ = [R{f}.happiness];
      idx = find(~isnan(happ));
      happ = zscore(happ(idx));
      for i = 1:length(idx)
        R{f}(idx(i)).happiness = happ(i);
      end
    end
  end % files
  
  % combine results structures?
  clear initialISIandSetup
  if size(files,1) > 1
    combineOrNot = 'combine';
    switch combineOrNot
      case 'combine' % then decide how many files to group
        if strcmpi(dataset,'Behav'), N = 3; else; N = 2; end
        for g = 1:size(files,1) / N % number of groups
          idx = [1:N] + (g-1)*N;
          RR{g} = horzcat(R{idx});
          subjNumbTemp(g) = subjNumbers(idx(1));
          % if need to create fMRI regressors, load the initial ISI and setup time
          if strcmpi(createfMRIreg,'Yes')
            idx = find(initialISIandSetupValues(:,1,1) == subjNumbTemp(g));
            if isempty(idx)
              warning(sprintf('no ISI values found for subject %d',subjNumbTemp(g)));
              initialISIandSetup(g,:) = NaN;
            else
              initialISIandSetup(g,:) = initialISIandSetupValues(idx,2:end); % go to the correct run and pick its number
            end
          end
        end
        R = RR;
        subjNumbers = subjNumbTemp;
      case 'don''t combine' % nothing to do
    end
  else; error('No data files found')
  end
  
  %% If Study 2, I can choose to analyse data from all the people with fMRI data (44 subjects), or all the people with good fMRI data (40 subjects)
  % We scanned 44 people, and excluded 4 that moved too much. I will report 
  % the behavioural data from 44 people, and the fMRI data from 40 people.
  if strcmpi(dataset,'fMRI')
    personsWithFMRIdata = [ 301 302 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 347 334 335 336 337 338 339 340 341 342 343 344 345 346 347 348 349]; 
    % personsWithGoodfMRIdata = [ 301 302 304 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320 321 322 323 324 325 326 327 328 329 330 331 332 333 347 334 335 336 337 338 339 340 ];
    subjNumbersOrig = subjNumbers;
    [subjNumbers,idx] = intersect(subjNumbers,personsWithFMRIdata);
    R = R(idx);
  end
  
  %% display number of subjects included:
  fprintf(1,' --- N subjects in this dataset: %d --- \n',length(subjNumbers))
  
  %% ------------------------------------------------------------------------
  %
  %     GOT DATA, NOW COLLECT AND DISPLAY HAPPINESS DATA - ANALYSIS l.1165+
  %
  % -------------------------------------------------------------------------
  
  %% Get happiness ratings for simple analysis - LMM analysis further below
  clear happiness* hap_idx* decisions decisionsCounts

  for i = 1:length(R) % for each subject
    happinessTemp = [R{i}(:).happiness];
    if strcmpi(zscoreSepEachRun,'No')
      idx = find(~isnan(happinessTemp));
      happ = zscore(happinessTemp(idx));
      happinessTemp(idx) = happ;
    end
    
    % analysis of the behavioural data
    % get the happiness values for the different outcomes
    hap_whoWon_soc = [1 1;1 0;0 1;0 0]; % first subject, then partner
    hap_whoWon_nonsoc = [1;0]; % only subject
    subjectWon = [R{i}(:).subjectWon];
    partnerWon = [R{i}(:).partnerWon];
    conditions = [R{i}(:).condition];
    happiness = [R{i}(:).happiness];
    dhappiness = NaN(size(happiness));  dhappiness(~isnan(happiness)) = [0 diff(happiness(2:2:end))];
    
    % Use difference in happiness for analysis
    if useDiffHappinessForSimplePlot
      happiness = dhappiness;
    end

    hap_idx = {};
    for ii = 1:length(hap_typeNames_soc) % for all 4 outcomes
      for cc = 1:2 % social and partner conditions
        hap_idx{ii,cc} = find(subjectWon == hap_whoWon_soc(ii,1) & partnerWon == hap_whoWon_soc(ii,2) & conditions == cc);
        if strcmpi(useAllTrialsForNonComp,'yes')
          try happiness_soc(i,ii,cc) = nanmean([happiness(hap_idx{ii,cc}) happiness(hap_idx{ii,cc}-1)]); catch happiness_soc(i,ii,cc) = NaN; end
        else
          try happiness_soc(i,ii,cc) = nanmean([happiness(hap_idx{ii,cc})]); catch happiness_soc(i,ii,cc) = NaN; end
        end
      end
    end
    for ii = 1:length(hap_typeNames_nonsoc)
      hap_idx{ii,1} = find(subjectWon == hap_whoWon_nonsoc(ii) & conditions==3);
      if strcmpi(useAllTrialsForNonComp,'yes')
        try happiness_nonsoc(i,ii,1) = nanmean([happiness(hap_idx{ii,1}) happiness(hap_idx{ii,1}-1)]); catch happiness_nonsoc(i,ii,1) = NaN; end
      else
        try happiness_nonsoc(i,ii,1) = nanmean([happiness(hap_idx{ii,1})]); catch happiness_nonsoc(i,ii,1) = NaN; end
      end
    end

  end % for each subject
  
  %% Report happiness ratings
  if sequentialFiguresForTalk; figsToDo = 1:6; else; figsToDo = 6; end
  for ff = figsToDo % for the different stages of the figure
    happiness = NaN(length(R),4,3); % dimensions: subjects; outcome [high|high, high|low, low|high, low|low]; condition [solo, social, partner]
    happinessNonsoc = NaN(length(R),4,3); % same size, just for the controls
    switch ff
      case 1 % first stage: only background
      case 2 % second stage: only non-soc
        happinessNonsoc(:,2:3,1) = happiness_nonsoc;
      case 3 % only add win-win
        happinessNonsoc(:,2:3,1) = happiness_nonsoc;
        happiness(:,1,2:3) = happiness_soc(:,1,:);
      case 4 % add win-lose
        happinessNonsoc(:,2:3,1) = happiness_nonsoc;
        happiness(:,1:2,2:3) = happiness_soc(:,1:2,:);
      case 5 % add lose-win
        happinessNonsoc(:,2:3,1) = happiness_nonsoc;
        happiness(:,1:3,2:3) = happiness_soc(:,1:3,:);
      case 6 % add the remainder
        happinessNonsoc(:,2:3,1) = happiness_nonsoc;
        happiness(:,:,2:3) = happiness_soc;
    end

    if happyFigClassicOrDotPlotFlag == 1
      figure('name','Happiness','position',[324   378   746   420]);
      if i>1 % more than 1 subject
        [ll1] = meansembar(1:3,happiness,1,niceColorsWarm(4),14,.2,'CI'); drawnow
        hold on
        [ll2] = meansembar(1:3,happinessNonsoc,1,[0 0 0;1 1 1;.5 .5 .5;0 0 0],14,.2,'CI'); drawnow
      else
        bar(squeeze(happiness)'); drawnow
      end
      if ~sequentialFiguresForTalk; legend([ll1 ll2(2:3)],[hap_typeNames_soc hap_typeNames_nonsoc],'location','eastoutside'); end
      ylabel('Happiness (z-scored)')
      set(gca,'position',[0.1300    0.2167    0.6998    0.7083],'ylim',[-1 1],'xlim',[.5 3.5],'xtick',1:3,'xticklabel',conditionNames([3 1 2]))
      grid on; box off; drawnow
    end
  end
  
  % For each subject, measure the impact of agency on happiness after negative outcomes for the partner
  responsibilityEffect = nanmean(squeeze(happiness(:,[2 4],2)-happiness(:,[2 4],3)),2)
  guiltEffect = responsibilityEffect; subjNumb = subjNumbers'; Tg = table(subjNumb,guiltEffect);
  writetable(Tg,[mfilepath 'csv/' dataset ' - BehavGuiltEffect.csv']);
  
  % Print happiness values to paper figure
  if happyFigClassicOrDotPlotFlag == 2
    if paperFigure
      figure(paperFig2)
      ah1 = axes('position',happyLotteryPanelPos+studyOffset);
    else
      figure('name',[dataset ' Happiness ratings dotPlots'],'position',[322   427   552   295]); ah1 = gca;
    end
    
    % Show happiness figure with dotplots
    % as soon as partner loses, participants feel better if partner took the decision
    axes(ah1)
    hh = []; S = []; happiness2 = [];
    happiness2(:,:,1:4) = permute(happiness(:,[1 3 2 4],2:3),[1 3 2]);
    S.colors = niceColorsWarm2(3); S.colors = S.colors(2:3,:);
    S.errorType = 'CI';
    S.markersize = 4;
    S.y = happiness2;
    hh = dotPlot(S);
    set(gca,'ylim',[-3.5 3.5],'ytick',-3:1:3,'xlim',[.2 8.8],'xtick',[1.5:2:8],'xticklabel',{'High\newlineHigh','Low\newlineHigh','High\newlineLow','Low\newlineLow',})
    box on
    hold on

    % add stars for "guilt effect": contrasting low partner outcomes in Social and Partner
    plot([5 6],[-3 -3],'k','linewidth',2)
    plot([7 8],[-3 -3],'k','linewidth',2)
    if strcmpi(dataset,'Behav') % p values are < 0.001 and = 0.002
      t_hdle=text(5.5,-2.8,'***'); t_hdle.HorizontalAlignment='center'; t_hdle.FontSize=20;
      t_hdle=text(7.5,-2.8,'**'); t_hdle.HorizontalAlignment='center'; t_hdle.FontSize=20;
    elseif strcmpi(dataset,'fMRI') % p values are = 0.01 and < 0.001
      t_hdle=text(5.5,-2.8,'**'); t_hdle.HorizontalAlignment='center'; t_hdle.FontSize=20;
      t_hdle=text(7.5,-2.8,'***'); t_hdle.HorizontalAlignment='center'; t_hdle.FontSize=20;
    end

    legend([hh(1:2).mu],{'Social','Partner'});
    ylabel('\Delta Happiness [std]')
    xlabel('Outcome participant (top) and partner (bottom)')
    title([studyName ': Happiness after lottery choices'])
    
  end
  
  %% ------------------------------------------------------------------------
  %
  %             NOW COLLECT, DISPLAY AND ANALYSE DECISIONS
  %
  % -------------------------------------------------------------------------
  %% Analyse decisions, but only for conditions 1 & 3, as an algorithm took decisions in condition 2
  if strcmpi(analyseDecisions,'Yes'); progressbar(0,0,'Decisions')
    EVdiff_str = {'EV_r_i_s_k_y - V_s_a_f_e'};
    Ntodo = length(R);
    for i = 1:Ntodo % for each subject
      conditions = [R{i}(:).condition];
      trialTypes = {'gain trials','mixed trials','loss trials'};
      for c = [1 3] % don't analyse computer decisions
        
        % ------ 0) Calculate basic variables ------
        idx_cond = find(conditions == c);
        Vsafe = [R{i}(idx_cond).safeOption]';
        temp = reshape([R{i}(idx_cond).riskyOption],2,[]);
        riskyHigh = temp(1,:)'; riskyLow = temp(2,:)';
        EVrisky = (riskyHigh + riskyLow) / 2;
        varRisky = riskyHigh - riskyLow;
        EVriskyMinusVsafe = EVrisky-Vsafe;
        chooseRisky = [R{i}(idx_cond).play]';
        choseSafe = 1-chooseRisky;
        
        % we sort the trials into gain, mixed, and loss trials, because we assume different risk preferences for them
        idx_gain = find(riskyLow >= 0 & Vsafe >= 0);
        idx_mixed = find(Vsafe == 0);
        idx_loss = find(riskyHigh <= 0 & Vsafe <= 0);
        trialType = [];  trialType(idx_gain) = 1;  trialType(idx_mixed) = 2;  trialType(idx_loss) = 3; trialType = trialType';
        
        trialType_indices = {idx_gain,idx_mixed,idx_loss};
        
        % Make a table with the variables we need
        Tdec = table(trialType, riskyHigh, riskyLow, EVrisky, varRisky, Vsafe, EVriskyMinusVsafe, chooseRisky, choseSafe);
        Tdec_allS{i,c} = Tdec; % gather all the decision data tables
        
        % ------ 1) Fit decisions as function of SVrisky-Vsafe with cumulative Gaussian ------
        if showIndivmodelNameFitFigs
          figure('name',sprintf('Subject %d',subjNumbers(i)));
        end
        
        if flag.fitBinnedEVdiffs
          [~,binEdges,ee] = histcounts(EVriskyMinusVsafe,binEdges);
          for iii = 1:length(binEdges)-1
            decisions(i,iii,c) = nanmean(chooseRisky(ee==iii));
            decisionsCounts(i,iii,c) = length(chooseRisky(ee==iii));
          end
          binCenters = (binEdges(1:end-1)+binEdges(2:end))/2;
          
          if showIndivmodelNameFitFigs
            if length(vars)>1; subplot(2,2,v); end
            hold on
            fitPsychFunction_js(binCenters/10,decisions(i,:,c),'plotFlag',gca,'functType',functType,'lowerBounds',lowerBounds,'upperBounds',upperBounds,'showText',0,'plotDottedLines','none','color',colors(c,:),'linewidth',2);
            title(EVdiff_str)
            %       ylabel('P(risky), \bar{x}\pm SEM','interpreter','tex')
            ylabel('P(risky) [M+-SEM]')
            grid on
          end
        end
        
        % ------ 2) Fit CARA model ------
        % now, new approach to model the choices people made.
        % We're using the constant absolute risk aversion (CARA) model to
        % obtain certainty equivalents for the risky options.
        
        if flag.findCARArhoInSearchSpace
          plotFlag = 0;
          for t = 1:length(trialTypes)
            idx = trialType_indices{t};
            CARA_rho(i,t,c) = findCARArhoInSearchSpace(riskyHigh(idx),riskyLow(idx),Vsafe(idx),chooseRisky(idx),rhoToTest,plotFlag);
          end
          CARA_rho_searchSpace(i,c) = findCARArhoInSearchSpace(riskyHigh,riskyLow,Vsafe,chooseRisky,rhoToTest,plotFlag);
        end
        
        % ------ 3) Fit probit function to model choices as function of difference in expected values between safe and risky:
        if flag.fitProbitFunction
          xForFit = -80:80; % delta EV (EVrisky - Vsafe) to use for display and calculation of risk premium
          yfit = {};
          [B,DEV,STATStemp] = glmfit(Tdec.EVriskyMinusVsafe,Tdec.chooseRisky,'binomial','link','probit');
          % Deviance test to see if this model fits better than one with just a constant term (see see https://ch.mathworks.com/help/stats/glmfit.html):
          [~,dev_noconstant] = glmfit(ones(size(Tdec.EVriskyMinusVsafe)),Tdec.chooseRisky,'binomial','link','probit','Constant','off');
          D = dev_noconstant - DEV; D = dev_noconstant - DEV; % D has a chi-square distribution with degrees of freedom equal to the difference in the number of estimated parameters in the model corresponding to DEV and the number of estimated parameters in the constant model.
          STATStemp.pModel = 1 - chi2cdf(D,length(B)-1); % Find the p-value for a deviance test.
          STATS_probitFitAllT(i,c) = STATStemp;
          yfit{1} = glmval(B,xForFit, 'probit');
          % Find 0-crossing point, which is a measure of the risk premium:
          [~,iii] = min((yfit{1}-.5).^2);
          riskPremiumProbitFit(i,c) = xForFit(iii);
        end
        drawnow
        
        % ------ 5) Calculate the SV of the risky option for each trial based
        % on CARA_rho, to see if choices are related to that now. Relate SV to happiness?
        if flag.findCARArhoInSearchSpace
          SVrisky = [];
          for t = 1:length(trialTypes)
            idx = Tdec.trialType == t;
            SVrisky(idx) = EU_CARA(Tdec.riskyHigh(idx),Tdec.riskyLow(idx),0.5,CARA_rho(i,t,c)); % subjective value of risky option under CARA given estimated rho
          end; SVrisky = SVrisky';
          Tdec.SVrisky = SVrisky;
          Tdec_allS{i,c}.SVrisky = SVrisky;
          SVdiff = Tdec.SVrisky-Tdec.Vsafe;
          Tdec.SVdiff = SVdiff;
          Tdec_allS{i,c}.SVdiff = SVdiff;
        end % flag.findCARArhoInSearchSpace
        
        % ------ 6) Fit decisions as function of SVrisky-SVsafe with cumulative Gaussian ------
        if flag.fitBinnedSVdiffs && ( flag.findCARArhoInSearchSpace )
          [~,binCenters,ee] = histcounts(SVdiff,binEdges);
          for iii = 1:length(binEdges)-1
            decisions(i,iii,c) = nanmean(chooseRisky(ee==iii));
          end
          binCenters = (binCenters(1:end-1)+binCenters(2:end))/2;
        end % flag.fitBinnedSVdiffs
        
      end % each experimental condition
      
      if showIndivmodelNameFitFigs; legend(conditionNames); end
      progressbar(i/Ntodo,0,'Decisions')
    end % each subject
  end % analyseDecisions
  
  %% Do participants choose on average the option with the higher EV?
  propRiskyChoicesIfEVdiffPos = [];
  for i = 1:Ntodo
    % first, define GLM design matrix
    X = [...
      Tdec_allS{i,1}.EVriskyMinusVsafe;...
      Tdec_allS{i,3}.EVriskyMinusVsafe...
      ];
    y = [Tdec_allS{i,1}.chooseRisky; Tdec_allS{i,3}.chooseRisky]; % risky choices
    propRiskyChoicesIfEVdiffPos(i) = mean([ y(X(:,1)>0); 1-y(X(:,1)<0) ]);
  end
  str = ttest_nice(propRiskyChoicesIfEVdiffPos-0.5);
  fprintf('Choices of option with higher EV: mean = %.1f, SD = %.1f, %s\n',mean(propRiskyChoicesIfEVdiffPos)*100,std(propRiskyChoicesIfEVdiffPos)*100,str)
  
  %% Check variations of CARA rho as a function of trialtype and expt. condition
  if flag.findCARArhoInSearchSpace
    disp('2-way ANOVA on CARA rho:')
    anova2_repmeas_autodesign(CARA_rho(:,:,[1 3]));
    disp('T-test on CARA rho obtained over all trials through search in space, comparing the 2 conditions:')
    ttest_nice(CARA_rho_searchSpace(:,1),CARA_rho_searchSpace(:,3));
  end
  
  %% Compare rho for self-for-both vs rho for self only
  figure('name',[dataset ' \rho across conditions (left=avg over trialtypes)'],'position',[221   240   815   345])
  rho = squeeze(mean(CARA_rho,2));
  subplot(1,2,1);
  regress_display(rho(:,1),rho(:,3),'regressType','robust ','inputNames',{'rho social','rho solo'},'axesHandle',gca);
  subplot(1,2,2);
  regress_display(makerow(CARA_rho(:,:,1)),makerow(CARA_rho(:,:,3)),'regressType','robust','inputNames',{'rho social','rho solo'},'axesHandle',gca);
  
  %% Fit generalized linear mixed-effects model to risky choices
  % Explanations here:
  % https://ch.mathworks.com/help/stats/generalized-linear-mixed-effects-models.html
  % https://ch.mathworks.com/help/stats/fit-a-generalized-linear-mixed-effects-model.html
  
  % 1) Make big table with data from all subjects
  Tdec = [];
  for i = 1:Ntodo
    temp1               = Tdec_allS{i,1}; % the self-for-both condition data
    temp1.subject       = i * ones(size(temp1.chooseRisky));
    temp1.condition     = ones(size(temp1.chooseRisky));
    temp2               = Tdec_allS{i,3}; % the self-only condition data
    temp2.subject       = i * ones(size(temp2.chooseRisky));
    temp2.condition     = 2 * ones(size(temp2.chooseRisky));
    Tdec = [Tdec; temp1; temp2];
  end
  % Mean-correct all terms first, to help interpret the coefficients.
  Tdec.condition      = Tdec.condition-1;
  Tdec.EVdiffMC       = ( Tdec.EVriskyMinusVsafe-mean(Tdec.EVriskyMinusVsafe) ); % !! do NOT change Euro cent into Euros
  Tdec.condMC         = Tdec.condition-mean(Tdec.condition);
  Tdec.EVcond_int_MC  = Tdec.condMC .* Tdec.EVdiffMC - mean(Tdec.condMC .* Tdec.EVdiffMC);
  
  % Write table out to verify model in R
  writetable(Tdec,[mfilepath 'csv/' dataset ' - Choices_singleTrialData.csv']);

  % 2) Make a simple linear mixed model with EVriskyMinusVsafe, condition and
  % their interaction as independent variables, subject as random term
  % factor, with random intercept.
  % NOTE June 2023: models with random slopes fit even better!
  % model = 'chooseRisky ~ 1 + EVdiffMC + condMC + EVcond_int_MC + (1 | subject)';
  % modelRS1 = 'chooseRisky ~ 1 + EVdiffMC + condMC + EVcond_int_MC + (1 + EVdiffMC | subject)';
  % modelRS2 = 'chooseRisky ~ 1 + EVdiffMC + condMC + EVcond_int_MC + (1 + EVdiffMC + condMC | subject)';
  model = 'chooseRisky ~ 1 + EVdiffMC * condMC + (1 + EVdiffMC + condMC | subject)';
  
  % 3) Fit a probit model
  disp('Fitting probit and linear probability models to Choice data...')
  glme = fitglme(Tdec,...
    model,...
    'Distribution','Binomial','Link','probit','FitMethod','Laplace',...
    'DummyVarCoding','effects');
  
  % 4) Fit a linear probability model, so that the coefficients can be better
  % interpreted
  lm = fitglme(Tdec,...
    model);
  
  % collect glme and lm for econ-style table
  choiceLMMs{currentStudy,1} = glme;
  choiceLMMs{currentStudy,2} = lm;

  %% 5) Make a figure with fitted data for the probit model
  %
  %   MUPRED = PREDICT(GLME) computes a vector MUPRED of predictions of the
  %   conditional mean of the response given the random effects at the
  %   original predictors used to create GLME. MUPRED includes contributions
  %   from both estimated fixed effects and approximate empirical Bayes
  %   predictors (EBPs) of the random effects.
  %
  %   MUPRED = PREDICT(GLME,T) uses predictor variables from the table T to
  %   predict the conditional mean of the response. T must contain all of the
  %   predictor variables used to create GLME.
  %
  %   [MUPRED,MUCI] = PREDICT(...) also returns the two-column matrix MUCI
  %   containing 95% pointwise confidence intervals for the predicted values.
  %   The lower limits of the bounds are in column 1, and the upper limits
  %   are in column 2.
  
  % Compute histograms of decisions at each EVdiff bin, for solo and social conditions
  clear hdle
  if paperFigure
    figure(paperFig1), axes('position',decPanelPos+studyOffset)
%     [hh,hhhh] = meansemplotWithDots(binCenters/100,decisions(:,:,[1 3])); hh(1).Visible = 'off'; hh(2).Visible = 'off'; hold on
%       colors = [0 0 .5; 0 0 0; .5 0 0]; % dark blue for solo and dark red for social
      colors = [.6 .6 1; 0 0 0; 1 .6 .6]; % dark blue for solo and dark red for social
      markers = '^ov'; % upper triangle for solo and lower triangle for social
    for c = [1 3]; for i=1:length(binCenters), [n,v]=hist(decisions(:,i,c)); keep=find(n); 
        scatter(binCenters(i)*ones(length(keep),1),v(keep),n(keep)*10,colors(c,:),'filled',markers(c));
        hold on, end, end
  else
    figure('Name',[dataset ' GLME fit'],'position',[ 360   427   346   271]);
  end
  
  % Print probit-model fits to figure
  % Nsubj = length(subjNumbers);
  Nd = 100; % N datapoints to model
  EVdiffMC = linspace(-50,50,Nd)';
  subject = zeros(Nd,1);
  colors = [0 0 1;1 0 0];
  
  for c = 1:2 % both conditions
    condMC = -1*(c-1.5) * ones(Nd,1);
    Tnew = table(EVdiffMC,condMC,subject);
    x = EVdiffMC;
    [mupred,muci] = predict(glme,Tnew); hold on
    hdle(c,:) = errorarea(x,mupred,muci,colors(c,:));
  end
  xlabel([EVdiff_str{1} '[Euro cents]'])
  ylabel('{\it p} (risky)')
  set(gca,'xlim',[-50 50])
  title([studyName ': Choices'])
  legend(hdle(:,1),conditionNames([3 1]),'location','NorthWest');
  
  %% If I assessed risk aversion using both CARA and the 0-crossing of the fitted probit function, I can check whether the derived measures are related:
  if flag.checkRiskAversionConsist && flag.fitProbitFunction && flag.findCARArhoInSearchSpace
    % check correlation between horizontal deviation in the probit fit
    % and rho
    %   riskPremiumProbitFit(abs(riskPremiumProbitFit)>50)=NaN;
    riskPrem = riskPremiumProbitFit / 100; % !! Transform Euro cent into Euro
    inputNames = {'\rho','Risk premium'};
    %   regress_display(CARA_rho(:),riskPremiumProbitFit_clean(:))
    figure('name',[dataset ' Consistency risk aversion measurements'],'position',[221   240   815   345])
    cc = [1 3]; % conditions to report on
    for c = 1:length(cc)
      subplot(1,2,c);
      riskPrem_temp = riskPrem(:,cc(c));
      rho_temp = mean(CARA_rho(:,:,cc(c)),2);
      regress_display(rho_temp,riskPrem_temp,...
        'axesHandle',gca,'inputNames',inputNames,'titlePrefix',conditionNames{cc(c)},'regressType','robust'); grid on
    end
    
    % Compare risk preference measure rho for Solo and Social:
    if paperFigure
      figure(paperFig1)
    else, 
      figure('name',[dataset ' Risk preference across conditions'],'position',[360   399   546   299]); 
    end
    
    facecolor = [.5 .5 .5]; withmdn = 0; style = 2;
    
    % Print estimated risk premium values to paper figure
    if paperFigure, figure(paperFig1), axes('position',riskPremPanelPos+studyOffset), else, subplot(1,2,2); end
    riskPrem_temp = riskPrem; riskPrem_temp(abs(riskPrem_temp)>.7) = 0;
    %   hh = distributionPlot(riskPrem_temp(:,[3 1]),'color',[.8 .8 .8],'showMM',4); hold on;
    %     hh{2}(1).LineWidth=2; hh{2}(2).LineWidth=2; hh{2}(1).Marker='.';
    violin(1,riskPrem_temp(:,3),'facecolor',facecolor,'withmdn',withmdn,'style',style,'scaling',1); hold on
    violin(2,riskPrem_temp(:,1),'facecolor',facecolor,'withmdn',withmdn,'style',style,'scaling',1); hold on
    %   hh{2}(1).Marker='o'; hh{2}(1).MarkerFaceColor=[1 0 0];
    %   hh{2}(2).Marker='^'; hh{2}(2).MarkerFaceColor=[0 1 0];
    plot(riskPrem_temp(:,[3 1])','k'); set(gca,'xtick',[1 2],'xticklabel',conditionNames([3 1]),'xlim',[.5 2.5]);
    m = mean(riskPrem_temp(:,[3 1]));
    CI = ci_t(riskPrem_temp(:,[3 1]),95); % find 95% confidence intervals for the mean, assuming normally-distributed data
    errL = m - CI(1,:);
    errH = CI(2,:) - m;
    % previously: err = std(riskPrem_temp(:,[3 1]))/sqrt(size(riskPrem_temp,1))
    errorbar([1 2],m,errL,errH,'.-r','linewidth',2);
    title([studyName ': Risk premium'])%,'position',[1    0.0027473])
    ylabel([EVdiff_str{1} ' [Euro]'])
    
    ttest_nice(riskPrem_temp(:,1),riskPrem_temp(:,3),'paired',1);

    % Print estimated risk aversion parameter rho values to paper figure
    if paperFigure, figure(paperFig1), axes('position',rhoPanelPos+studyOffset), else, 
      subplot(1,2,1); 
    end
    violin(1,rho(:,3),'facecolor',facecolor,'withmdn',withmdn,'style',style,'scaling',1); hold on
    violin(2,rho(:,1),'facecolor',facecolor,'withmdn',withmdn,'style',style,'scaling',1); hold on
    plot(rho(:,[3 1])','k'); set(gca,'xtick',[1 2],'xticklabel',conditionNames([3 1]),'xlim',[.5 2.5])
    m = mean(rho(:,[3 1]));
    CI = ci_t(rho(:,[3 1]),95); % find 95% confidence intervals for the mean, assuming normally-distributed data
    errL = m - CI(1,:);
    errH = CI(2,:) - m;
    % previously: errorbar([1 2],mean(rho(:,[3 1])),std(rho(:,[3 1]))/sqrt(size(rho,1)),'.-r','linewidth',2);
    errorbar([1 2],m,errL,errH,'.-r','linewidth',2);
    title([studyName ': \rho'],'fontsize',24)
    ylabel('risk seeking \leftrightarrow risk averse')
    
    ttest_nice(rho(:,1),rho(:,3),'paired',1);
    
  end
  
  %% ------------------------------------------------------------------------
  %
  %        ANALYSED DECISIONS, NOW COMPUTATIONAL MODELING OF HAPPINESS
  %
  % -------------------------------------------------------------------------
  
  %% Rutledge-like computational analysis
  for i = 1:length(R) % for each subject
    % first, prepare the variables for the computational analysis
    choseSafe = [R{i}(:).play];
    %   choseGamble = 1-choseSafe;
    try, nonsocialTrial = ([R{i}(:).condition] == 3);
    catch; nonsocialTrial = ~isnan([R{i}(:).amountThisTrialPartner]); end
    try, socialActive = ([R{i}(:).condition] == 1);
    catch; socialActive = [R{i}.condition]; end
    try, socialPassive = ([R{i}(:).condition] == 2); end
    amountReceivedBySubject = [R{i}.amountThisTrialSubject];
    CR = amountReceivedBySubject.*choseSafe;
    EV = mean(reshape([R{i}(:).riskyOption],2,[])) .* (1-choseSafe);
    RPE = EV-amountReceivedBySubject.*(1-choseSafe);
    
    happinessTemp = [R{i}(:).happiness];
    if strcmpi(zscoreSepEachRun,'No')
      idx = find(~isnan(happinessTemp));
      happ = zscore(happinessTemp(idx));
      happinessTemp(idx) = happ;
    end
    
    matrix = zeros(length(R{i}),20);
    matrix(:,1) = [R{i}(:).condition];
    matrix(:,3) = [R{i}(:).safeOption];
    matrix(:,4:5) = reshape([R{i}(:).riskyOption],2,[])';
    matrix(:,7) = [R{i}(:).play];
    matrix(:,8) = [R{i}(:).amountThisTrialSubject];
    matrix(:,10) = happinessTemp;
    matrix(:,13) = isnan([R{i}(:).amountThisTrialPartner])+1;
    matrix(:,16) = [R{i}(:).amountThisTrialPartner];
    matrix(:,17) = [R{i}(:).Agl_Abs];
    matrix(:,18)= [R{i}(:).MinEV];
    matrix(:,19)= [R{i}(:).ChoiceEV];
    matrix(:,20)= [R{i}(:).Ch_diff];
    matrix(:,21) = [R{i}(:).PE];
    
    switch whichTrials
      case 'all'
        idx = 1:size(matrix,1);
      case 'social'
        idx = find(~nonsocialTrial);
      case 'non-social'
        idx = find(nonsocialTrial);
      case 'social-active'
        idx = find(socialActive);
      case 'social-passive'
        idx = find(socialPassive);
      case 'social-passive and non-social'
        idx = [find(socialPassive) find(nonsocialTrial)];
      case 'social-active and non-social'
        idx = [find(socialActive) find(nonsocialTrial)];
    end
    matrix = matrix(idx,:);
    
    alldata(i).socdata{1} = matrix;
    %   alldata(i).Ntrials = [R{i}.Ntrials];
  end % subject
  
  
  %% Run Rutledge-like computational analysis
  clear compResults mtxs ff
  for whatToDo = compToDoList % 1 = fit the data, 2 = create fMRI regressors; 0 = parameter recovery; done after having fitted real data.
    disp('Fitting computational models for participant number: ')
    % first, fit behavioural data, then create fMRI regressors if desired
    
    for n=1:length(alldata), % n subjects
      fprintf(1,'%d,',n);
      for m=1:length(alldata(n).socdata), % m runs?
        % temp is the data, rows are trials, columns are variables:
        % 3: safe values
        % 4 to 5: risky values
        % 7: 1 if gambled, 0 if chose safe option
        % 8: reward obtained by self
        % 10: happiness score
        % 13: 1 if participant chose, 2 if other person chose
        % 16: other person reward
        % 17: Agl_Abs
        % 18: MinEV
        % 19: ChoiceEV
        % 20: Ch_diff
        % 21: PE
        eval(sprintf('temp = alldata(%d).socdata{%d};',n,m)); temp(1,10:12)=nan; %collect matrices and toss first rating
        t2 = temp(:,10); rawhappy = t2(~isnan(t2)); happyind = find(~isnan(t2)); temp(isnan(temp(:,16)),16)=0; %other reward=0 when solo trial
        nrate=length(rawhappy); ntrial=size(temp,1); x=zeros(nrate,ntrial);
        evmtx=x; certainmtx=x; rpemtx=x; rewardmtx=x; othercertainmtx=x; otherrewardmtx=x; otherrpemtx=x; ipemtx=x; envymtx=x; guiltmtx=x; notenvymtx=x; notguiltmtx=x; youchosemtx=x; theychosemtx=x;
        Agl_Absmtx=x; MinEVmtx=x; ChoiceEVmtx=x; Ch_diffmtx=x; PEmtx=x;
        switch whatToDo
          case 0 % parameter recovery; same trials as for fitting
            trialsToFit = happyind;
          case 1 % fit the data
            trialsToFit = happyind;
          case 2 % create fMRI regressors
            trialsToFit = 1:ntrial;
        end
        
        for q = 1:length(trialsToFit)
          temp2            = temp(1:trialsToFit(q),:); %clip out all trials up to rating
          tempev           = mean(temp2(:,4:5),2) .* double(temp2(:,7)==1); %0 if no gamble or error, ev if gambled
          temprpe          = temp2(:,8) .* double(temp2(:,7)==1) - tempev; %0 if no gamble or error, rpe if gambled
          tempreward       = temp2(:,8) .* double(temp2(:,7)==1); % gamble rewards (0 if chose certain)
          tempcertain      = temp2(:,3) .* double(temp2(:,7)==0); % 0 if gambled, otherwise certain amount
          tempothercertain = temp2(:,16) .* double(temp2(:,7)==0) .* double(temp2(:,13)==2); % other reward if no gamble
          tempotherreward  = temp2(:,16) .* double(temp2(:,7)==1) .* double(temp2(:,13)==2); % other reward if gamble
          tempotherrpe     = temp2(:,16) .* double(temp2(:,7)==1) - tempev; % rpe for other
          tempipe          = (sign(temprpe))   .* double(temp2(:,1)<3) .* (tempotherreward-tempev); %sign or your rpe x sign of theirs
          tempenvy         = (sign(temprpe)<0) .* double(temp2(:,1)<3) .* (tempotherreward-tempreward); %you lost and they got this much more than you
          tempguilt        = (sign(temprpe)>0) .* double(temp2(:,1)<3) .* (tempreward-tempotherreward); %you won and you got this much more than them
          tempnotenvy      = (sign(temprpe)<0) .* double(temp2(:,1)<3) .* double(~(tempotherreward-tempreward)); %you lost and they also lost dummy
          tempnotguilt     = (sign(temprpe)>0) .* double(temp2(:,1)<3) .* double(~(tempreward-tempotherreward)); %you won and they also won dummy
          tempyouchose     = double(temp2(:,13)==1);
          temptheychose    = double(temp2(:,13)==2);
          evmtx(q,1:length(tempev))                     = fliplr(transpose(tempev)) / 100; %convert to ?
          certainmtx(q,1:length(tempcertain))           = fliplr(transpose(tempcertain)) / 100;
          rpemtx(q,1:length(temprpe))                   = fliplr(transpose(temprpe)) / 100;
          rewardmtx(q,1:length(tempreward))             = fliplr(transpose(tempreward)) / 100;
          ipemtx(q,1:length(tempipe))                   = fliplr(transpose(tempipe)) / 100;
          othercertainmtx(q,1:length(tempothercertain)) = fliplr(transpose(tempothercertain)) / 100;
          otherrewardmtx(q,1:length(tempotherreward))   = fliplr(transpose(tempotherreward)) / 100;
          otherrpemtx(q,1:length(temprpe))              = fliplr(transpose(tempotherrpe)) / 100;
          envymtx(q,1:length(tempenvy))                 = fliplr(transpose(tempenvy)) / 100;
          guiltmtx(q,1:length(tempguilt))               = fliplr(transpose(tempguilt)) / 100;
          notenvymtx(q,1:length(tempnotenvy))           = fliplr(transpose(tempnotenvy)); %dummy
          notguiltmtx(q,1:length(tempnotguilt))         = fliplr(transpose(tempnotguilt)); %dummy
          youchosemtx(q,1:length(tempyouchose))         = fliplr(transpose(tempyouchose)); %1 if you chose
          theychosemtx(q,1:length(temptheychose))       = fliplr(transpose(temptheychose)); %1 if they chose
          Agl_Abs=double(temp2(:,17));
          MinEV=double(temp2(:,18));
          ChoiceEV=double(temp2(:,19));
          Ch_diff=double(temp2(:,20));
          PE=double(temp2(:,21));
          Agl_Absmtx(q,1:length(Agl_Abs))               = fliplr(transpose(Agl_Abs)) / 100; %convert to ?
          MinEVmtx(q,1:length(MinEV))                   = fliplr(transpose(MinEV)) / 100; %convert to ?
          ChoiceEVmtx(q,1:length(ChoiceEV))             = fliplr(transpose(ChoiceEV)) / 100; %convert to ?
          Ch_diffmtx(q,1:length(Ch_diff))               = fliplr(transpose(Ch_diff)) / 100; %convert to ?
          PEmtx(q,1:length(PE))                         = fliplr(transpose(PE)) / 100; %convert to ?
        end;
        clear mtx; mtx.certainmtx=single(certainmtx); mtx.evmtx=single(evmtx); mtx.rpemtx=single(rpemtx); mtx.rewardmtx=single(rewardmtx);
        mtx.ipemtx=single(ipemtx); mtx.othercertainmtx=single(othercertainmtx); mtx.otherrewardmtx=single(otherrewardmtx); mtx.otherrpemtx=single(otherrpemtx);
        mtx.envymtx=single(envymtx); mtx.guiltmtx=single(guiltmtx); mtx.youchosemtx=single(youchosemtx); mtx.theychosemtx=single(theychosemtx);
        mtx.notenvymtx=single(notenvymtx); mtx.notguiltmtx=single(notguiltmtx);
        mtx.Agl_Abs=single(Agl_Absmtx); mtx.MinEV=single(MinEVmtx); mtx.PE=single(PEmtx);
        mtx.ChoiceEV=single(ChoiceEVmtx); mtx.Ch_diff=single(Ch_diffmtx);
        mtx.conds=temp(1:trialsToFit(end),1);
        eval(sprintf('alldata(%d).mtx{%d}=mtx;',n,m));
        eval(sprintf('alldata(%d).rawhappy{%d}=rawhappy;',n,m));
        eval(sprintf('temp = alldata(%d).socdata{%d};',n,m)); temp(1,10:12)=nan; %toss first rating
        t2 = temp(:,10); zhappy = t2(~isnan(t2));
        eval(sprintf('alldata(%d).zhappy{%d}=zhappy;',n,m))
        
        % do computational modelling
        if whatToDo == 1 %2 & strcmpi(createfMRIreg,'Yes')
          if strcmpi(compModelName,'ActPass')
            compResults(n,1) = fit_happy_model_inequalityActPass(mtx,zhappy,NaN);
          elseif strcmpi(compModelName,'otherRPEactPass')
            compResults(n,1) = fit_happy_model_otherRPEactPass(mtx,zhappy,NaN);
          elseif strcmpi(compModelName,'GuiltEnvy')
            compResults(n,1) = fit_happy_model_guiltenvy(mtx,zhappy,NaN);
          elseif strcmpi(compModelName,'InequalSimple')
            compResults(n,1) = fit_happy_model_inequalitySimple(mtx,zhappy,NaN);
          elseif strcmpi(compModelName,'Basic')
            compResults(n,1) = fit_happy_model_basic(mtx,zhappy,NaN);
          elseif strcmpi(compModelName,'Responsibility')
            compResults(n,1) = fit_happy_model_responsibility(mtx,zhappy,NaN);
          elseif strcmpi(compModelName,'ResponsibilityRedux')
            compResults(n,1) = fit_happy_model_responsibility_redux(mtx,zhappy,NaN);
          elseif strcmpi(compModelName,'GuiltEnvy+UtilityFunction')
            compResults(n,1) = fit_happy_model_guiltenvy_utl(mtx,zhappy,NaN);
          elseif strcmpi(compModelName,'Compare')
            compResults(n,1) = fit_happy_model_basic(mtx,zhappy,NaN);
            compResults(n,2) = fit_happy_model_inequalitySimple(mtx,zhappy,NaN);
            compResults(n,3) = fit_happy_model_guiltenvy(mtx,zhappy,NaN);
            compResults(n,4) = fit_happy_model_responsibility(mtx,zhappy,NaN);
            compResults(n,5) = fit_happy_model_responsibility_redux(mtx,zhappy,NaN);
          end
          if showIndivmodelNameFitFigs
            for m = 1:size(compResults,2)
              figure;
              subplot(2,1,1)
              plot(zhappy,'b');
              hold on; plot(compResults(n,m).happypred,'r')
              dataNames = {'happiness','predicted happiness'};
              legend(dataNames);
              title([R{n}(1).subjName ' ' whichTrials], 'interpreter','none')
              ylabel(sprintf('R2 = %.2f',compResults(n,m).r2));
              subplot(2,1,2)
              bar(compResults(n,m).b)
              xlabel_oblique(compResults(n,m).paramNames)
              title([compResults(n,m).compModelName ' - parameters'], 'interpreter','none')
              drawnow
            end
          end
          %     regress_display(zhappy,results(n).happypred,'inputNames',dataNames)
        end % whatToDo == 1
        
        % use tau to create fMRI regressors if desired
        if whatToDo == 2 && strcmpi(createfMRIreg,'Yes')
          idx = strmatch('decay',compResults(n,1).paramNames); % look for the decay parameter
          tau = compResults(n,1).b(idx); %decay constant
          if useTauForFMRIreg
            [fMRIregressors,fMRIregressor_names] = create_fMRIregressors_otherRPEactPass(tau,mtx,R{n},initialISIandSetup(n,:)); % use results of the fit to create fMRI regressors
            %           [fMRIregressors,fMRIregressor_names] = create_fMRIregressors(tau,mtx,R{n},initialISIandSetup(n,:)); % use results of the fit to create fMRI regressors
          else
            [fMRIregressors,fMRIregressor_names] = create_fMRIregressors(0,mtx,R{n},initialISIandSetup(n,:)); % use results of the fit to create fMRI regressors
          end
          q=R{n}(1).subjName; us = findstr(q,'_');
          subjNo = q(us(5)+1:us(6)-1);
          if double(subjNo(1)) > 57
            subjNo = q(us(4)+1:us(5)-1);
          end
          for ses = 1:length(fMRIregressors)
            dlmwrite(['fMRIregressors_' subjNo '_' num2str(ses) '.txt'],fMRIregressors{ses})
            figure; imagesc(fMRIregressors{ses});
            set(gca,'xtick',[1:size(fMRIregressors{ses},2)]+.4)
            xlabel_oblique(fMRIregressor_names(1:size(fMRIregressors{ses},2)),15,12,1);
            set(gca,'position',[0.1300    0.1800    0.7750    0.7450])
            ylabel('Trial from start of run1')
            title(['fMRIregressors for subject ' subjNo ', run ' num2str(ses)])
          end % fMRI ses
        end % fMRI regressors
      end % runs
      mtxs{n} = mtx;
      zhappyAll(n).zhappy = zhappy;
      subjectNos(n).subjectNo = n*ones(length(zhappy),1);
    end % subjects
  end % whatToDo
  fprintf(1,'done\n');
  
  % Make table with computational results, write it out
  if compModelling && whatToDo ~= 0 % so no parameter recovery
    for m = 1:size(compResults,2)
      predictedHappiness  = vertcat(compResults(:,m).happypred);
      realHappiness       = vertcat(zhappyAll(:).zhappy);
      subjectNo           = vertcat(subjectNos(:).subjectNo);
      ThapComp            = table(realHappiness,predictedHappiness,subjectNo);
      writetable(ThapComp,[mfilepath 'csv/' dataset ' - Happiness and happiness predicted by ' compResults(1,m).modelName '.csv'])
    end
  end
  if ~compModelling || whatToDo == 0 % then read already computed predicted happiness values
    ThapComp = readtable([mfilepath 'csv/' dataset ' - Happiness and happiness predicted by ResponsibilityRedux.csv']);
  end
  % Show fit quality of model
  if strcmpi(study,'both (-> paper)')
    figure(paperFig2)
    axes('position',happyPredPanelPos+studyOffset);
  else
    figure
  end
  regress_display(ThapComp.realHappiness,ThapComp.predictedHappiness,'axesHandle',gca,'inputNames',{'Reported \Delta happiness [std]','Pred. \Delta happiness [std]'},'regressType','robust')
  title([studyName ': Fit of comp. model'])
  
  %% Show results of computational modelling
  try disp('Individual fits of the models tested (R2):')
    disp(reshape([compResults.r2],n,[])); end
  
  if exist('compResults','var')
    
    % combine the results in a matrix
    clear resMtx
    f1 = figure('name','RL model parameters','position',[440   139   560   659]);
    for m = 1:size(compResults,2)
      for n = 1:size(compResults,1)
        resMtx{m}(n,:) = [subjNumbers(n) compResults(n,m).r2 compResults(n,m).r2adj compResults(n,m).aic compResults(n,m).bic compResults(n,m).b]; % with aic and bic
      end
      resMtxVarnames{m} = [{'Subject number', 'R2', 'R2adj', 'aic', 'bic'}, compResults(n,m).paramNames];
      
      % display
      if length(alldata) > 1
        % group results
        figure(f1); subplot(size(compResults,2),1,m)
        showWhat = [6:size(resMtx{m},2)];
        plot([.5 length(showWhat)+.5],[0 0],'k')
        try violinDots(resMtx{m}(:,showWhat));
        catch, meansemplot(resMtx{m}(:,showWhat));
        end
        try xlabel_oblique(resMtxVarnames{m}(showWhat)); end
        title(compResults(1,m).modelName)
        
        % individual variation
        figure('name',['Individual variation of parameters for compModelName ' compResults(1,m).modelName]);
        try
          subplot(4,1,1)
          bar(resMtx{m}(:,2)); title('R2')
          subplot(4,1,2)
          bar(responsibilityEffect); title('Responsibility effect: Effect of self - other''s decision on happiness in other loss trials')
          subplot(4,1,3)
          bar(resMtx{m}(:,end-1)); title(compResults(n,m).paramNames(end-1))
          subplot(4,1,4)
          bar(resMtx{m}(:,end)); title(compResults(n,m).paramNames(end))
        catch
          subplot(4,1,1)
          plot(resMtx{m}(:,2)); title('R2')
          subplot(4,1,2)
          plot(responsibilityEffect); title('Responsibility effect: Effect of self - other''s decision on happiness in other loss trials')
          subplot(4,1,3)
          plot(resMtx{m}(:,end-1)); title(compResults(n,m).paramNames(end-1))
          subplot(4,1,4)
          plot(resMtx{m}(:,end)); title(compResults(n,m).paramNames(end))
        end
        % write out data to analyse in SPSS or such
        dlmwrite(['ComputFitData_' compResults(n,m).modelName '.csv'],resMtx{m})
      end
      disp('Mean and median R2:')
      disp([mean(resMtx{m}(:,2)) median(resMtx{m}(:,2))])
    end
  else; disp('No computational modelling done')
  end
  
  if strcmpi(compModelName,'compare') % then make a table showing the comparison results and do some stats on the parameters
    modelNames = {compResults(1,:).modelName};
    Nmodels = length(modelNames);
    figure('name','RL model comparison','position',[954   480   446   325])
    for m = 1:Nmodels
      aic(m) = sum(resMtx{m}(:,3));
      bic(m) = sum(resMtx{m}(:,4));
    end
    subplot(2,1,1); bar(bic); title('BIC'); set(gca,'xtick',1:Nmodels)
    subplot(2,1,2); bar(aic); title('AIC'); set(gca,'xtick',1:Nmodels,'xticklabel',modelNames)
    
    % do statistics on the model parameters:
    disp('Testing if the model parameters CR, EV, RPE and decay are above 0, and if EV < RPE:')
    for m = 1:Nmodels; disp(modelNames{m});
      for tt = 5:8; ttest_nice(resMtx{m}(:,tt)); end;
      ttest_nice(resMtx{m}(:,6),resMtx{m}(:,7));
    end

    % save the calculated model parameters for each model, and do stats on the Responsibilty model parameters:
    for m = 1:Nmodels
      params = vertcat(compResults(:,m).b);
      paramNames = compResults(1,m).paramNames;
      ParamFitTable = array2table(params,'VariableNames',paramNames);
      writetable(ParamFitTable,[mfilepath 'csv/' dataset ' - ' modelNames{m} ' - fittedParameters.csv']);
      
      % stats on the Responsibiltiy model parameters:
      if m == 4 % the Responsibility model
        disp('Comparing the model parameters sRPE > social_pRPE:')
        str = signrank_nice(ParamFitTable.sRPE, ParamFitTable.social_pRPE); disp(str)

        disp('Comparing the model parameters social_pRPE > 0:')
        str = signrank_nice(ParamFitTable.social_pRPE); disp(str)

        disp('Comparing the model parameters social_pRPE > partner_pRPE:')
        str = signrank_nice(ParamFitTable.social_pRPE, ParamFitTable.partner_pRPE); disp(str)

      elseif m == 5 % the ResponsibilityRedux model
        disp('Comparing the model parameters sRPE > social_pRPE:')
        str = signrank_nice(ParamFitTable.sRPE, ParamFitTable.social_pRPE); disp(str)

        disp('Comparing the model parameters social_pRPE > 0:')
        str = signrank_nice(ParamFitTable.social_pRPE); disp(str)

      end
    end

    % make paper table
    varIdxToReport = [2 3 5 4]; % R2, adjR2, BIC, AIC
    for m = 1:length(modelNames)
      dataForTable(m,:) = [mean(resMtx{m}(:,2)) mean(resMtx{m}(:,3)) round(sum(resMtx{m}(:,5))) round(sum(resMtx{m}(:,4)))];
    end
    Tcomp = array2table(dataForTable,'VariableNames', resMtxVarnames{1}(varIdxToReport));
    Tcomp.modelName = modelNames';
    writetable(Tcomp,[mfilepath 'csv/' dataset ' - computational models of happiness.xlsx']);

    % do some stats to compare the model fits
    disp('Results of likelihood ratio test comparing the different models:')
    % run likelihood ratio tests (LR = -2 * ln(L0/L1) = -2 * (logL1 - logL0)), follows Chi-squared distribution, with DF = diff in N params:
    % See for example here: https://www.geeksforgeeks.org/r-language/likelihood-ratio-test/
    comparisons = [[1 4];[2 4];[3 4];[5 4]]; % the numbers of the models to be compared: compare Responsibility (model 4) to all others
    for pair = 1:size(comparisons,1)
      m0 = comparisons(pair,1); m1 = comparisons(pair,2); % m0 is the simpler, m1 the more complex model
      % LR = -2 * ([compResults(:,m0).logL] - [compResults(:,m1).logL]); % for each participant separately
      LR = -2 * (sum([compResults(:,m0).logL]) - sum([compResults(:,m1).logL])); % summing over participants
      df = length(compResults(1,m1).b) - length(compResults(1,m0).b);
      pval = chi2test(LR,df); % LRs follow a Chi-square distribution
      fprintf('Is model "%s" better than model "%s"? LikelihoodRatioTest = %.2f, p = %f\n', modelNames{m1}, modelNames{m0}, LR, pval)
    end

  end % compare models
  
  if strcmpi(compModelName,'ParamRecovery')
    %% Run parameter recovery analysis, only for ResponsibilityRedux model. Need to load real parameters and compute mtx matrices (done above).
    
    % First, load the real fitted parameters
    paramDir = 'csv/';
    params = readtable([paramDir dataset ' - ResponsibilityRedux - fittedParameters.csv']);
    
    % Next, for each subject of each study, simulate data based on real
    % parameters, fit those data, and check what parameters come out. Vary the
    % noise (units: SD of happiness). Display results.
    noiseLevel = 1; Nrep = 10; displayFlag = 0;
    if ~(length(mtxs) == size(params,1)); error('Data and parameter matrices have different sizes - check which study we''re processing here'); end
    progress_report(0,'Parameter recovery')
    recov_params_all = [];
    for sub = 1:length(mtxs)
      progress_report(sub/length(mtxs),'Parameter recovery')
      mtx = mtxs{sub};
      paramsThisSub = table2array(params(sub,:));
      recov_params = parameterRecovery_happy_model_responsibility_redux(mtx,paramsThisSub,noiseLevel,Nrep,displayFlag);
      recov_params_all(:,:,sub) = recov_params; % dimensions: Nrep, paramNo, subNo.
    end % for each subject
    save(['bin/compModels/Recovered model parameters for ' dataset ' study ' num2str(Nrep) 'iter.mat'],'recov_params_all','params')

    % Now, display
    figure('Name',['Model parameter recovery for ' dataset ' study']);
    for p = 1:size(params,2)
      subplot(2,3,p);
      % plot(table2array(params(:,p)),squeeze(mean(recov_params_all(:,p,:))),'.') % x axis: real parameters across subjects, y axis: recovered parameters
      plot(repmat(table2array(params(:,p)),1,Nrep)',squeeze(recov_params_all(:,p,:)),'.') % x axis: real parameters across subjects, y axis: recovered parameters
      title(params.Properties.VariableNames{p},'interpreter','none')
    end
  end % parameter recovery

  %% ------------------------------------------------------------------------
  %
  %        NOW LINEAR MIXED MODEL ANALYSIS OF HAPPINESS DATA
  %
  % -------------------------------------------------------------------------
  data = [];
  for i = 1:length(R) % for each subject
    for j = 1:length(R{i}) % for each trial
      dat = [...
        R{i}(j).happiness...
        R{i}(j).condition...
        R{i}(j).safeOption/100 ...
        R{i}(j).riskyOption/100 ...
        R{i}(j).play...
        R{i}(j).subjectWon...
        R{i}(j).partnerWon...
        R{i}(j).amountThisTrialSubject/100 ...
        R{i}(j).amountThisTrialPartner/100 ...
        i...
        ];
      data = [data; dat];
    end
  end
  Thap = array2table(data,'VariableNames',{'happiness' 'cond' 'Vsafe' 'riskyHi' 'riskyLo' 'choseRisky' 'subjectWon' 'partnerWon' 'rewardSubj' 'rewardPart' 'subject'});
  Thap.EV           = ( ( Thap.riskyHi + Thap.riskyLo ) / 2 ) .* Thap.choseRisky + ( Thap.Vsafe .* ( 1 - Thap.choseRisky) );
  ineq              = NaN(size(data,1),1); ineq(Thap.rewardSubj ~= Thap.rewardPart) = 1; ineq(Thap.rewardSubj == Thap.rewardPart) = 0; ineq(isnan(Thap.rewardPart)) = NaN; % social trials in which participants got different outcomes vs same outcomes
  ineqAdv           = NaN(size(data,1),1); ineqAdv(Thap.subjectWon==1 & Thap.partnerWon==0) = 1;    ineqAdv(Thap.subjectWon==1 & Thap.partnerWon==1) = 0;     % trials in which participant only won vs. both won
  ineqDisadv        = NaN(size(data,1),1); ineqDisadv(Thap.subjectWon==0 & Thap.partnerWon==1) = 1; ineqDisadv(Thap.subjectWon==1 & Thap.partnerWon==1) = 0;  % trials in which partner only won vs. both won
  Thap.ineq         = ineq;
  Thap.ineqAdv      = ineqAdv;
  Thap.ineqDisadv   = ineqDisadv;
  Thap.sRPE         = Thap.rewardSubj - Thap.EV;
  Thap.pRPE         = Thap.rewardPart - Thap.EV;
  Thap.subjDecided  = Thap.cond ~= 2;
  Thap.socialTrial  = Thap.cond <  3;
  
  % Write table out to verify model in R - remove trials with happiness =
  % NaN, trials were the safe option was chosen, and non-social trials
  Thap2 = Thap(~isnan(Thap.happiness) & ~isnan(Thap.partnerWon),:);
  writetable(Thap2,[mfilepath 'csv/' dataset ' - Happiness_singleTrialData_socialRiskyChoicesOnly.csv']);
  
  % Make a nice simple plot relating reward to happiness change
  if strcmpi(study,'both (-> paper)')
    figure(paperFig2)
    axes('position',happyRewardPanelPos1+studyOffset); regress_display(Thap.rewardSubj, Thap.happiness,'axesHandle',gca,'colors',colors(1,:),'inputNames',{'Participant reward [Euro]','\Delta Happiness [std]'},'regressType','robust');
    title({[studyName ':'],'Impact of participant','reward on happiness'})
    axes('position',happyRewardPanelPos2+studyOffset); regress_display(Thap.rewardPart, Thap.happiness,'axesHandle',gca,'colors',colors(2,:),'inputNames',{'Partner reward [Euro]','\Delta Happiness [std]'},'regressType','robust');
    % title('   Impact of partner\newlinereward on happiness')
    title({[studyName ':'],'Impact of partner','reward on happiness'})
  else
    figure
    subplot(1,2,1); regress_display(Thap.rewardSubj, Thap.happiness,'axesHandle',gca,'colors',colors(1,:),'inputNames',{'Participant reward [Euro]','\Delta Happiness [std]'},'regressType','robust');
    title('Impact of participant\newlinereward on happiness')
    subplot(1,2,2); regress_display(Thap.rewardPart, Thap.happiness,'axesHandle',gca,'colors',colors(2,:),'inputNames',{'Partner reward [Euro]','\Delta Happiness [std]'},'regressType','robust');
    title('   Impact of partner\newlinereward on happiness')
  end
  
  % Testing simple models relating happiness to reward
  model = {}; modelName = {}; lm_data = [];
  model{1} = 'happiness ~ 1 + rewardSubj + rewardPart + (1|subject)';         % modelName{1} = 'Subject and partner reward';
  model{2} = 'happiness ~ 1 + rewardSubj + rewardPart + ineq + (1|subject)';  % modelName{2} = 'Subject reward + inequality (lotteries only)';
  lm = {};
  for m = 1:length(model)
    lm{m} = fitglme(Thap,model{m});
  end
  for m = 2:length(model)
    results = compare(lm{1},lm{m})%,'CheckNesting',true)
  end
  [~,best] = min(results.AIC);
  disp('========= Linear mixed model regressing happiness on reward; best model was:')
  lm{best}
  printEconomicsStyleRegressionTableForGLMs(lm, modelName, [dataset ' - LMM happiness on reward'], [mfilepath 'csv/' dataset ' - LMM happiness on reward.csv']);
  
  % show residuals of best-fitting model
  [~,best] = min(results.AIC);
  figure('position',[319   137   889   798],'name','Residuals of LM on happiness, simple and responsibility models');
  [H,p,jbstat] = jbtest(lm{best}.residuals);
  if H; str = 'Residuals not normally distributed'; else; str = 'Residuals normally distributed'; end
  tit = sprintf('Jarque-Bera stat = %.1f, p = %.3f',jbstat,p);
  subplot(2,2,1); plot(Thap.happiness,lm{best}.residuals,'.'); axis tight; xlabel('Happiness'); ylabel('Residuals'); title(str)
  subplot(2,2,2); qqplot(lm{best}.residuals); xlabel('Normal quantiles'); ylabel('Residuals'); title(tit)

  % Condition models
  model = {}; lm = {}; modelName = {};
  model{1} = 'happiness ~ 1 + cond + (1|subject)';
  model{2} = 'happiness ~ 1 + subjDecided + (1|subject)';
  model{3} = 'happiness ~ 1 + socialTrial + (1|subject)';
  model{4} = 'happiness ~ 1 + subjDecided + socialTrial + (1|subject)';
  for m = 1:length(model)
    lm{m} = fitglme(Thap,model{m});
    disp([lm{m}.LogLikelihood lm{m}.Rsquared.Adjusted]);
  end
  for m = 2:length(model)
    results = compare(lm{1},lm{m})%,'CheckNesting',true)
  end
  [~,best] = min(results.AIC);
  disp('========= Linear mixed model regressing happiness on subjDecided and social or not; best model was:')
  lm{best}
  printEconomicsStyleRegressionTableForGLMs(lm, modelName, [dataset ' - LMM happiness on conditions'], [mfilepath 'csv/' dataset ' - LMM happiness on conditions.csv']);

  % Responsibility models
  model = {}; lm = {}; modelName = {};
  model{1} = 'happiness ~ 1 + subjectWon + partnerWon + subjDecided + (1|subject)';
  model{2} = 'happiness ~ 1 + subjectWon + partnerWon + subjectWon:partnerWon + subjDecided + subjDecided:subjectWon + subjDecided:partnerWon + (1|subject)';
  for m = 1:length(model)
    lm{m} = fitglme(Thap,model{m});
    disp([lm{m}.LogLikelihood lm{m}.Rsquared.Adjusted]);
  end
  for m = 2:length(model)
    results = compare(lm{1},lm{m})%,'CheckNesting',true)
  end
  disp('========= Linear mixed model regressing happiness after lottery outcomes based on who won and who decided; best model was:')
  lm{end}

  % show residuals of best-fitting Responsibility LMM model
  [~,best] = min(results.AIC);
  [H,p,jbstat] = jbtest(lm{best}.residuals);
  if H; str = 'Residuals not normally distributed'; else; str = 'Residuals normally distributed'; end
  tit = sprintf('Jarque-Bera stat = %.1f, p = %.3f',jbstat,p);
  subplot(2,2,3); plot(Thap.happiness,lm{best}.residuals,'.'); axis tight; xlabel('Happiness'); ylabel('Residuals'); title(str)
  subplot(2,2,4); qqplot(lm{best}.residuals); xlabel('Normal quantiles'); ylabel('Residuals'); title(tit)

  % show effect of outcome for partner X decision-maker on happiness
  q = lm{end}.Coefficients(7,:);
  fprintf('Study %d, happiness, linear model t statistic for %s: t(%d) = %.2f, p value = %f, beta = %.3f, 95%% CI = [%.3f %.3f]\n', currentStudy, q{1,1}, q{1,5}, q{1,4}, q{1,6}, q{1,2}, q{1,7}, q{1,8})

  printEconomicsStyleRegressionTableForGLMs(lm, {}, [dataset ' - LMM happiness on outcomes and responsibility'], [mfilepath 'csv/' dataset ' - LMM happiness on outcomes and responsibility.csv']);
  
  %% ------------------------------------------------------------------------
  %
  %                   NOW INTER-INDIVIDUAL DIFFERENCES 
  %         EXPLORATION, NOT REPORTED IN PAPER BECAUSE TOO LOW N
  %
  % -------------------------------------------------------------------------
  
  %% Relate responsibilityEffect to some outcomes of the computational analysis
  if strcmpi(modelName,'Compare')
    q=[responsibilityEffect resMtx{4}(:,9:10) resMtx{4}(:,9)-resMtx{4}(:,10)];
    idx=~isnan(responsibilityEffect);
    names={'responsibility','otherRPEact','otherRPEpass','diff'};
    disp(names)
    [rr,pp]=corr(q(idx,:))
  end
  
  %% Check if measures of risk aversion relate to the responsibility effect:
  if flag.fitProbitFunction && flag.findCARArhoInSearchSpace
    figure('name','Does risk aversion relate to the responsibility effect?','position',[202    83   840   612])
    respo = 'Respons. effect happiness';
    subplot(2,2,1);
    inputNames = {respo,'\rho'};
    regress_display(responsibilityEffect,mean(CARA_rho(:,:,1),2),'regressType','robust','axesHandle',gca,'inputNames',inputNames)
    subplot(2,2,2);
    inputNames = {respo,'Risk premium'};
    regress_display(responsibilityEffect,mean(riskPremiumProbitFit(:,1),2),'regressType','robust','axesHandle',gca,'inputNames',inputNames)
    subplot(2,2,3);
    inputNames = {respo,'\Delta \rho (social-solo)'};
    regress_display(responsibilityEffect,mean(CARA_rho(:,:,1)-CARA_rho(:,:,3),2),'regressType','robust','axesHandle',gca,'inputNames',inputNames)
    subplot(2,2,4);
    inputNames = {respo,'\Delta Risk premium (social-solo)'};
    regress_display(responsibilityEffect,mean(riskPremiumProbitFit(:,1)-riskPremiumProbitFit(:,3),2),'regressType','robust','axesHandle',gca,'inputNames',inputNames)
  end
  
end % which dataset

%% If output for paper, report mixed-model regressions on choices and make econ-style mixed model output tables

if paperFigure
  disp('Stats for probit models of choices')
  for st = 1:2 % Study 1 or 2

    % Stats for effect of EVdiff in probit model of choices
    q = choiceLMMs{st,1}.Coefficients(2,:); % Probit model
    fprintf('Study %d, probit model t statistic for %s: t(%d) = %.2f, p value = %d, beta = %.3f, 95%% CI = [%.3f %.3f]\n', st, q{1,1}, q{1,5}, q{1,4}, q{1,6}, q{1,2}, q{1,7}, q{1,8})

    % Stats for effect of condition (social vs solo)
    q = choiceLMMs{st,1}.Coefficients(3,:); % Probit model for Study 1
    fprintf('Study %d, probit model t statistic for %s: t(%d) = %.2f, p value = %d, beta = %.3f, 95%% CI = [%.3f %.3f]\n', st, q{1,1}, q{1,5}, q{1,4}, q{1,6}, q{1,2}, q{1,7}, q{1,8})
  end

  % Print all model stats to file
  modelNames = {'Study 1 probit','Study 1 linear','Study 2 probit','Study 2 linear'};
  printEconomicsStyleRegressionTableForGLMs(choiceLMMs([1 3 2 4]), modelNames, 'Choices', [mfilepath 'csv/LMMs choices.csv']);
end

%% If paperFigure, add letters and print as SVG
panLetterFontSize = 30; % 40 previously
if paperFigure
  figure(paperFig1)
  axes('position',[0 0 1 1],'visible','off')
  text(repmat([.02 .49 .73],1,2),kron([.89 .44],[1 1 1]),cellstr(['ABCDEF']'),'fontsize',panLetterFontSize)
  cd([mfilepath 'Figures/'])
  printSVG
  printfig
  figure(paperFig2)
  axes('position',[0 0 1 1],'visible','off')  %.41
  text(repmat([.02 .23 .44 .67],1,2),kron([.89 .44],[1 1 1 1]),cellstr(['ABCDEFGH']'),'fontsize',panLetterFontSize)
  printSVG
  printfig
  cd(mfilepath)
end

%% Show parameter recovery
% A nice supplementary materials figure of the parameter recovery results
compModelParameterRecoveryFigure

%% POSSIBLE TO DOs
set(0,'defaultAxesFontSize',12) % back to normal defaults

disp('Could check whether the variance of the lotteries correlates with their EV or EV diff to alternative safe option, as both variance and deltaEV can affect choices')
disp('Could relate sympathy levels to responsibilityEffect, if I had more participants')
