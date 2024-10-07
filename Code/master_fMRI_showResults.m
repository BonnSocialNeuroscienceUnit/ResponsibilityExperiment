% SHOW FIGURES OF MAIN FMRI RESULTS AND CALCULATE EFFECT SIZES FROM BETAS
function master_fMRI_showResults
%
% js 2019-2024

% general settings
panLetterFontSize = 28;
sliceFontSize = 18;
labelFontSize = 16;
LMMcoeffBarColor = [.8 .8 .8];
tFlg = 0; % if 1, adds text title to panels
figName1 = 'Figure 4 - fMRIresults';
figName2 = 'Figure 5 - fMRI connectivity results'; % will be split into Fig 5 and Supp Fig 2 by hand after saving

% panel settings
i = 1; PP{i}.label = 'A'; PP{i}.name = 'Choose risky > safe';               PP{i}.position =  [0.012  0.68   0.24   0.28];
i = 2; PP{i}.label = 'B'; PP{i}.name = 'Decision social > solo';            PP{i}.position =  [0.27    0.68   0.42   0.28];
i = 3; PP{i}.label = 'C'; PP{i}.name = 'LMM decision ROIs';                 PP{i}.position1 =  [0.77   0.79   0.10   0.17];  PP{i}.position2 =  [0.89   0.79   0.10   0.17];
i = 4; PP{i}.label = 'D'; PP{i}.name = 'Outcome risky > safe';              PP{i}.position1 = [0.002  0.35   0.26   0.28];   PP{i}.position2 = [0.12   0.35   0.38  0.28];
i = 5; PP{i}.label = 'E'; PP{i}.name = 'LMM outcome ROIs';                  PP{i}.position1 =  [0.54   0.46   0.10   0.17];  PP{i}.position2 =  [0.66   0.46   0.10   0.17];
i = 6; PP{i}.label = 'F'; PP{i}.name = 'Guilt effect';                      PP{i}.position =  [0.75   0.35    0.28   0.28];
i = 7; PP{i}.label = 'G'; PP{i}.name = 'CR+EV TD model';                    PP{i}.position =  [0.002  0.02   0.26   0.28];
i = 8; PP{i}.label = 'H'; PP{i}.name = 'Other RPE Act>Pass TD model';       PP{i}.position = [0.27 0.02 0.3 0.28];
i = 9; PP{i}.label = 'I'; PP{i}.name = 'Model parameters in left STS';      PP{i}.position = [0.65 0.05 0.29 0.22];

% directories
mfilepath             = fileparts(mfilename('fullpath'));
addpath(genpath(mfilepath)) % add code and subfolders to path
baseDir               = [fileparts(mfilepath) '/fMRIresults/']
baseDirDecisions      = [baseDir 'decision/'];
baseDirOutcomes       = [baseDir 'outcome/'];
baseDirModelBased     = [baseDir 'model-based/'];
baseDirConnect        = [baseDir 'PPI/'];
figuresDir            = [mfilepath '/Figures'];

% reference anatomical scan
anat = [baseDir 'mni152_2009bet.nii.gz'];

% whole-brain results images
decisionRiskySafe         = [baseDirDecisions 'risky>safe_0p05FWE_clust.nii'];
decisionSocialSolo        = [baseDirDecisions 'social>solo_0p05FWE_clust.nii'];
outcomeRiskySafe          = [baseDirOutcomes 'risky>safe_0p05FWE_clust.nii'];
guiltEffect               = [baseDirOutcomes 'guiltEffect_0p05FWE_SVC_aIns.nii'];
modelBased_CRandEV        = [baseDirModelBased 'CR+EV_0p001u_k70_0p05FWE.nii'];
modelBased_pRPEactVsPass  = [baseDirModelBased 'pRPEsocial>pRPEpartner_0p001u_k70.nii'];
PPI_aInsDecision          = [baseDirConnect 'aIns_seed/Risky>SafeXSolo>Social_0p001u_k30_2IFGs.nii'];
PPI_STSdecision           = [baseDirConnect 'STS_seed/leftIFG_Risky>SafeXSocial>Solo_0p001u_k20.nii'];

% display
if isempty(which('spm')); spm12; end

%% Figure with main activation findings
figure('name',figName1,'position',[42  12  1263  696])
for p = 1:length(PP) % [1 2 4 6 7]
  % -------------------------------------------------------------------------
  %% A)
  if p == 1
    pan = PP{p};
    axes('position',pan.position)
    imagescSPM2imgsHighResMNI(decisionRiskySafe,anat,[],10,'coronal',sliceFontSize,gca)
    text(-55,1,pan.label,'fontsize',panLetterFontSize); if tFlg; title(pan.name); end
  end
  % Calculate Cohen's d for the 3 clusters while we're here
  disp('Cohen''s d risky vs safe decisions at the following locations:')
  clust = {'Group response at -10 8 -6 VStriaL.mat','Group response at 10 12 -4 VStriaR.mat'};
  for c = 1:length(clust)
    load([baseDirDecisions 'BetasDecisionROIsForLMM/' clust{c}])
    dec_risky = mean(RFXresponse.response(:,[[7 9] [7 9]+19]),2); % mean of Decision risky social and Decision risky non-social, both sessions
    dec_safe = mean(RFXresponse.response(:,[[4 6] [4 6]+19]),2); % mean of Decision safe social and Decision safe non-social, both sessions
    d = cohens_d(dec_risky-dec_safe); % this is a within-subject, one-sample contrast
    fprintf('%s: d = %.2f\n',clust{c},d)
  end

  % -------------------------------------------------------------------------
  %% B)
  if p == 2
    pan = PP{p};
    axes('position',pan.position)
    imagescSPM2imgsHighResMNI(decisionSocialSolo,anat,[],[-40 0],'sagittal',sliceFontSize,gca)
    text(-55,1,pan.label,'fontsize',panLetterFontSize); if tFlg; title(pan.name); end
  end
  % Calculate Cohen's d for the 3 clusters while we're here
  disp('Cohen''s d social vs solo decisions at the following locations:')
  clust = {'Group response at 0 -62 38 Precun.mat','Group response at -40 -54 26 TPJ_L.mat','Group response at 6 52 22 DMPFC.mat'};
  for c = 1:length(clust)
    load([baseDirDecisions 'BetasDecisionROIsForLMM/' clust{c}])
    dec_soc = mean(RFXresponse.response(:,[[4 7] [4 7]+19]),2); % mean of Decision safe social + Decision risky social, both sessions
    dec_nsoc = mean(RFXresponse.response(:,[[6 9] [6 9]+19]),2); % mean of Decision safe non-social + Decision risky non-social, both sessions
    d = cohens_d(dec_soc-dec_nsoc); % this is a within-subject, one-sample contrast
    fprintf('%s: d = %.2f\n',clust{c},d)
  end
  % -------------------------------------------------------------------------
  %% C)
  if p == 3
    pan = PP{p};
    load('variablesLMMdecisionForMainFigure.mat')
    roiNames = fieldnames(coefsAllROIs);
    idx = [4 5]; % the ROIs to show
    ll = 3; xpos = [1:ll ll+1.5 ll+3];

    axes('position',pan.position1);
    coefs = coefsAllROIs.Precuneus; coefCIs = coefCIsAllROIs.Precuneus;
    bar(xpos,coefs,'FaceColor',LMMcoeffBarColor);
    hold on; errorbar(xpos,coefs,coefCIs(1,:),coefCIs(2,:),'.k','linewidth',1);
    set(gca,'xtick',[1:ll ll+1.5 ll+3],'xlim',[.2 ll+3.8],'FontSize',labelFontSize)
    xlabel_oblique(coefNames)
    title(roiNames{idx(1)})
    ylabel({'LMM coefficients';'+- 95%CI'});
    text_relLoc(-.74,1,pan.label,'black',panLetterFontSize); if tFlg; title(pan.name); end
    grid on

    axes('position',pan.position2);
    coefs = coefsAllROIs.TPJ; coefCIs = coefCIsAllROIs.TPJ;
    bar(xpos,coefs,'FaceColor',LMMcoeffBarColor);
    hold on; errorbar(xpos,coefs,coefCIs(1,:),coefCIs(2,:),'.k','linewidth',1);
    set(gca,'xtick',[1:ll ll+1.5 ll+3],'xlim',[.2 ll+3.8],'FontSize',labelFontSize)
    xlabel_oblique(coefNames)
    title(roiNames{idx(2)})
    grid on
  end
  % -------------------------------------------------------------------------
  %% D)
  if p == 4
    pan = PP{p};
    axes('position',pan.position1);
    imagescSPM2imgsHighResMNI(outcomeRiskySafe,anat,[],4,'sagittal',sliceFontSize,gca,0)
    text(-55,1,pan.label,'fontsize',panLetterFontSize); if tFlg; title(pan.name); end
    axes('position',pan.position2);
    imagescSPM2imgsHighResMNI(outcomeRiskySafe,anat,[],[-3 44],'transversal',sliceFontSize,gca) ;% I want Z = -3
  end
  % Calculate Cohen's d for these clusters while we're here
  disp('Cohen''s d risky vs safe outcomes at the following locations:')
  clust = {'-30 20 -10 InsL','30 22 -10 InsR','2 42 36 DMPFC','48 -24 -8 MidTemp R','-10 0 -6 VStriaL',...
    '8 4 2 VStriaR','40 -48 46 ParietalR','42 38 24 FrontalR','-46 -44 46 ParietalL'};
  for c = 1:length(clust)
    load([baseDirOutcomes 'betasRiskySafe/Group response at ' clust{c} '.mat'])
    out_safe = mean(RFXresponse.response(:,[[10:12] [10:12]+19]),2); % mean of Outcome safe active, passive, nsoc, both sessions
    out_risky = mean(RFXresponse.response(:,[[13:18] [13:18]+19]),2); % mean of Outcome pos and neg, both sessions
    d = cohens_d(out_risky-out_safe); % this is a within-subject, one-sample contrast
    fprintf('Group response at %s: d = %.2f\n',clust{c},d)
  end
  % -------------------------------------------------------------------------
  %% E)
  if p == 5
    pan = PP{p};
    load('variablesLMMoutcomeForMainFigure.mat')
    roiNames = fieldnames(coefsAllROIs);
    idx = [3 4]; % the ROIs to show
    ll = 3; xpos = [1:ll ll+1.5 ll+3];

    axes('position',pan.position1);
    coefs = coefsAllROIs.InsulaL; coefCIs = coefCIsAllROIs.InsulaL;
    bar(xpos,coefs,'FaceColor',LMMcoeffBarColor);
    hold on; errorbar(xpos,coefs,coefCIs(1,:),coefCIs(2,:),'.k','linewidth',1);
    set(gca,'xtick',[1:ll ll+1.5 ll+3],'xlim',[.2 ll+3.8],'FontSize',labelFontSize)
    xlabel_oblique(coefNames)
    title(roiNames{idx(1)})
    ylabel({'LMM coefficients';'+- 95%CI'})
    text_relLoc(-.66,1,pan.label,'black',panLetterFontSize); if tFlg; title(pan.name); end
    grid on

    axes('position',pan.position2);
    coefs = coefsAllROIs.InsulaR; coefCIs = coefCIsAllROIs.InsulaR;
    bar(xpos,coefs,'FaceColor',LMMcoeffBarColor);
    hold on; errorbar(xpos,coefs,coefCIs(1,:),coefCIs(2,:),'.k','linewidth',1);
    set(gca,'xtick',[1:ll ll+1.5 ll+3],'xlim',[.2 ll+3.8],'FontSize',labelFontSize)
    xlabel_oblique(coefNames)
    title(roiNames{idx(2)})
    grid on
  end
  % -------------------------------------------------------------------------
  %% F)
  if p == 6
    pan = PP{p};
    axes('position',pan.position);
    imagescSPM2imgsHighResMNI(guiltEffect,anat,[],24,'coronal',sliceFontSize,gca,1)
    text(-55,1,pan.label,'fontsize',panLetterFontSize); if tFlg; title(pan.name); end

    % Calculate Cohen's d for this cluster while we're here
    disp('Cohen''s d guilt effect at the following location:')
    clust = {'-28 24 -4 LeftInsulaGuiltEffect'};
    load([baseDirOutcomes 'BetasGuiltEffect/Group response at ' clust{1} '.mat'])
    out_risky_neg_participant = mean(RFXresponse.response(:,[16 16+19]),2); % mean of Outcome risky neg social, both sessions
    out_risky_neg_partner = mean(RFXresponse.response(:,[17 17+19]),2); % mean of Outcome risky neg partner, both sessions
    d = cohens_d(out_risky_neg_participant-out_risky_neg_partner); % this is a within-subject, one-sample contrast
    fprintf('Group response at %s: d = %.2f\n',clust{1},d)
  end
  % -------------------------------------------------------------------------
  %% G)
  if p == 7
    pan = PP{p};
    axes('position',pan.position)
    imagescSPM2imgsHighResMNI(modelBased_CRandEV,anat,[],9,'coronal',sliceFontSize,gca)
    text(-55,1,pan.label,'fontsize',panLetterFontSize); if tFlg; title(pan.name); end
  end

  % Calculate Cohen's d for this cluster while we're here
  disp('Cohen''s d CR+EV at the following locations:')
  clust = {'-14 8 -8 VStriaL_CR+EV_0p001u','10 10 -4 VStriaR_CR+EV_0p001u'};
  for c = 1:length(clust)
    load([baseDirModelBased 'Group response at ' clust{c} '.mat'])
    CR_EV = mean(RFXresponse.response(:,[[2:3] [2:3]+6]),2); % mean of CR and EV, both sessions
    d = cohens_d(CR_EV); % this is a within-subject, one-sample contrast
    fprintf('Group response at %s: d = %.2f\n',clust{c},d)
  end
%% H)
  if p == 8
    pan = PP{p};
    axes('position',pan.position)
    imagescSPM2imgsHighResMNI(modelBased_pRPEactVsPass,anat,[],-52,'sagittal',sliceFontSize,gca,1)
    text(-55,1,pan.label,'fontsize',panLetterFontSize); if tFlg; title(pan.name); end
  end

  % Calculate Cohen's d for this cluster while we're here
  disp('Cohen''s d of greater partner RPE following participant than partner decisions at the following location:')
  clust = {'-52 -32 0 leftSTS_0p05FWEclus'};
  load([baseDirModelBased 'Group response at ' clust{1} '.mat'])
  pRPE_participant = mean(RFXresponse.response(:,[5 11]),2); % mean of pRPE participant decisions, both sessions
  pRPE_partner = mean(RFXresponse.response(:,[6 12]),2); % mean of pRPE partner decisions, both sessions
  d = cohens_d(pRPE_participant-pRPE_partner); % this is a within-subject, one-sample contrast
  fprintf('Group response at %s: d = %.2f\n',clust{1},d)
  
  %% I)
  if p == 9
    pan = PP{p};
    axes('position',pan.position)
    load([baseDirModelBased 'Group response at -52 -32 0 STSleft_pRPE_social>partner.mat']);
    temp = RFXresponse.response(:,[5:6 11:12]); % data = reshape(temp,80,2);
    data = []; data(:,:,1)=temp(:,1:2); data(:,:,2)=temp(:,3:4); % 3rd dimension are sessions, 2nd dimension are conditions
    xlab = {'pRPE participant','pRPE partner'}; % xlab = char(RFXresponse.condNames(5:6));
    [~,~,legHdle] = meansemplotDots(data,'errorType','CI','markersize',24); legHdle.String = {'Sess 1','Sess 2'};
    % violinDots(1:2,data,.8); % alternative display
    set(gca,'XTickLabel',xlab,'FontSize',labelFontSize)
    ylabel({'Model-predicted';'response (A.U.)'})
    grid on
    text_relLoc(-.24,1.14,pan.label,'black',panLetterFontSize); title(pan.name); 
  end
  
end

% Print figure to file
cd(figuresDir)
printSVG
printfig
cd(mfilepath)

%% figure with connectivity results at decision time, interaction RiskySafeXSocialSolo
figure('name',figName2,'position',[757   572   665   422])
i = 10; PP{i}.label = 'A'; PP{i}.name = 'Connectivity with left insula'; PP{i}.position1 =  [0.05  0.55   0.5   0.37]; PP{i}.position2 =  [0.6  0.59   0.3   0.33];
i = 11; PP{i}.label = 'B'; PP{i}.name = 'Connectivity with left STS';    PP{i}.position1 =  [0.05  0.08   0.5   0.37];  PP{i}.position2 =  [0.6  0.12   0.3   0.33]; 

% Insula seed: Show IFG cluster on slice of anatomical brain
pan = PP{10};
axes('position',pan.position1)
imagescSPM2imgsHighResMNI(PPI_aInsDecision,anat,[],46,'sagittal',sliceFontSize,gca,1)
text(-75,1,pan.label,'fontsize',panLetterFontSize); if tFlg; title(pan.name); end

% Insula seed: Load parameter estimates and display
load([baseDirConnect 'aIns_seed/Group response at 46 16 22 rightIFG_0p05FWEclust.mat'])
myresp = []; % collect responses and order in 2x2 design
myresp(:,1:2,1) = RFXresponse.response(:,[22 20]); % safe solo, safe social
myresp(:,1:2,2) = RFXresponse.response(:,[25 23]); % risky solo, risky social
axes('position',pan.position2);
[~,~,legHdle] = meansemplotDots(myresp,'markersize',24);
legHdle.String = {'Safe','Risky'}; legHdle.Location='southeast';
set(gca,'fontsize',labelFontSize,'xtick',[1 2],'xticklabel',{'Solo','Social'});
grid on

% Insula seed: Calculate Cohen's d for the right IFG cluster
fprintf('Cohen''s d of connectivity change between left insula seed and right IFG cluster as function of condition: ')
interaction = -myresp(:,1,1)+myresp(:,2,1)+myresp(:,1,2)-myresp(:,2,2); % -safe solo + safe social + risky solo - risky social
d = cohens_d(interaction); % this is a within-subject, one-sample contrast
fprintf('d = %.2f\n',d)

% STS seed: Show IFG cluster on slice of anatomical brain
pan = PP{11};
axes('position',pan.position1)
imagescSPM2imgsHighResMNI(PPI_STSdecision,anat,[],-48,'sagittal',sliceFontSize,gca,1)
text(-75,1,pan.label,'fontsize',panLetterFontSize); if tFlg; title(pan.name); end

% STS seed: Load parameter estimates and display
load([baseDirConnect 'STS_seed/Group response at -48 14 6 leftIFG_Risky>SafeXSocial>Solo.mat']);
myresp = []; % collect responses and order in 2x2 design
myresp(:,1:2,1) = RFXresponse.response(:,[22 20]); % safe solo, safe social
myresp(:,1:2,2) = RFXresponse.response(:,[25 23]); % risky solo, risky social
axes('position',pan.position2);
[~,~,legHdle] = meansemplotDots(myresp,'markersize',24);
legHdle.String = {'Safe','Risky'}; legHdle.Location='southeast';
set(gca,'fontsize',labelFontSize,'xtick',[1 2],'xticklabel',{'Solo','Social'});
grid on

% STS seed: Calculate Cohen's d for the left IFG cluster
fprintf('Cohen''s d of connectivity change between STS seed and left IFG cluster as function of condition: ')
interaction = myresp(:,1,1)-myresp(:,2,1)-myresp(:,1,2)+myresp(:,2,2); % safe solo - safe social - risky solo + risky social
d = cohens_d(interaction); % this is a within-subject, one-sample contrast
fprintf('d = %.2f\n',d)

% Print figure to file
cd(figuresDir)
printSVG
printfig
cd(mfilepath)
