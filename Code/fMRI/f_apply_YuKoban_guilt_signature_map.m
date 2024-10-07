% Script from Tor Wager's canlab github to compare an image to a "neural
% signature", a multivariate predictive pattern, for example for guilt
% ("Guilt-Related Brain Signature", or GRBS, from Yu, Koban et al., CerCor
% 2020) or reward ("A neural signature for reward", Speer et al.,
% Neuroimage 2023).

%% Set up directories and input data
baseDir = '/Users/johannesschultz/sciebo/CURRENT WORK/Responsibility/';

% the published guilt-associated activation pattern
% The guilt weight map is Yu_guilt_SVM_sxpo_sxpx_EmotionForwardmask.nii.gz in the repository.
guilt_pattern_Yu = [baseDir 'Code/bin/Yu_guilt_SVM_sxpo_sxpx_EmotionForwardmask.nii'];

% group response
my_responses = [baseDir 'fMRIresults/outcome/guiltEffectEachPartic.nii']; % 4D image with each subject's guilt response as separate map

% If you have our "CanlabCore" repository from Github on your Matlab path, 
% an easy way to apply the signature pattern to a new dataset is the "apply_mask" 
% method for fmri_data objects.  For example, this code will load a sample dataset 
% and apply the signature, and return a vector of one pattern response per input image.

%% Add CanLab code to path
canlabcode = {'/Volumes/LaCie_non-RAID/Canlab/','/Users/johannesschultz/Documents/MatlabToolboxes/CanlabCore-master'};
for c = 1:length(canlabcode)
  if exist(canlabcode{c},'dir'), addpath(genpath(canlabcode{c})); end
end

%% Load image to test
subjNumb = [ 101 102 103 104 106 110 120 124 125 126 129 136 155 166 171 173 174 175 44 55 62 63 65 69 70 71 74 77 79 81 83 84 85 90 92 93 94 95 96 98 ]'; % all subjects in the order stored in the 4D image
my_responses_fmri_data_object = fmri_data(my_responses, [], 'noverbose');  % loads images as fmri_data object

%% Load reference image to compare to
% guilt = load_image_set('guilt');  % Load the Guilt Behavior map as an fmri_data object
% myfile = which('Yu_guilt_SVM_sxpo_sxpx_EmotionForwardmask.nii.gz');
guilt_Yu_fmri_data_object = fmri_data(guilt_pattern_Yu, [], 'noverbose');  % load the Guilt Behavior map as an fmri_data object

%% Compare test image(s) to a reference image; returns the dot product between each test image and the reference
guiltResponseSimilarity = apply_mask(my_responses_fmri_data_object, guilt_Yu_fmri_data_object, 'pattern_expression'); % Apply the weight map, returns vector of dot products

% Put the results into a table
Tf = table(subjNumb,guiltResponseSimilarity);

fprintf('Is the neural guilt response significantly higher than 0? ')
if kstest(guiltResponseSimilarity) % results not normally distributed
  [str,p_guiltSignature] = signtest_nice(guiltResponseSimilarity); % tests if median is different from 0
else
  [str,~,p_guiltSignature] = ttest_nice(guiltResponseSimilarity);
end
fprintf([str '\n'])
% figure; plot(rewardResponse, 'o'); title('Reward response scores for each participant'); % that's a control

% Compute effect size: Area under ROC and Cliff's Delta, using measures-of-effect-size-toolbox (Hentschke and Stüttgen 2017)
addpath(genpath('/Users/johannesschultz/Documents/MatlabToolboxes/measures-of-effect-size-toolbox-master'))
guiltResponseEffectSize = mes(guiltResponseSimilarity,0,'auroc'); % Area under the ROC is a classic effect size measure for non-parametric data, see Documentation_MESToolbox.pdf page 20
auroc = guiltResponseEffectSize.auroc;
cliff = (2 * auroc)-1; % See Documentation_MESToolbox.pdf page 21: Cliff and Auroc are linearly related
fprintf('Effect size (AUROC and Cliff''s Delta) of the guilt response: %.2f, %.2f\n',auroc,cliff)

%% Correlate neural guilt activation pattern response to my behavioural guilt effect
% Load behavioural data, put into table
behavGuiltEffectFilename = [baseDir 'Behav data/BehavDuringfMRI/fMRI - GuiltOrResponsibilityEffect.csv'];
Tb = readtable(behavGuiltEffectFilename);

% Now match the subject numbers in the behavioural and fMRI data
[subjInBothDatasets,IB,IF] = intersect(Tb.subjNumb,Tf.subjNumb); % [C,IA,IB] = intersect(A,B) returns index vectors IA and IB such that C = A(IA) and C = B(IB).
behav = Tb.guiltEffect(IB);
brain = Tf.guiltResponseSimilarity(IF);

% Test correlation (Spearman's Rho)
[rho,pval] = corr(behav,brain,'Type','Spearman');
fprintf('Is the Behavioural guilt effect related to the GBRS dot product? No: rho=%.3f, p=%.3f\n',rho,pval)

%% Display
figure('name','Apply Yu et al 2020 Guilt Signature','position',[18 354 1324 379]);
subplot(1,5,[1 2]); plot(guiltResponseSimilarity,'o');
title('Dot-product with Yu et al 2020 GRBS');
set(gca,'xtick',1:length(subjNumb)) % ,'xticklabel',subjNumb) % don't show subject ID for anonymity
xlabel('Participant No.')

subplot(1,5,3)
violinDots(guiltResponseSimilarity);
title(sprintf('Median > 0 signtest: p = %.3f', p_guiltSignature))
set(gca,'XTickLabel','Similarity (dot product with GRBS)')
box on

subplot(1,5,[4 5])
regress_display(behav,brain,'inputNames',{'Behavioural guilt effect','GBRS dot product'},'axesHandle',gca);
