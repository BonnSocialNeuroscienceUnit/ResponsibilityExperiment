%% Makes a nice supplementary materials figure of the parameter recovery results
printFigFlag = 0; % set to 1 to print figure in 2 png and svg formats

% Load the results of the parameter recovery: parameters obtained by
% fitting synthetic data created with fitted parameters, with 1SD of noise
% added to the happiness data. 1000 iterations per subject.
% baseDir = '/Users/johannesschultz/sciebo/CURRENT WORK/Responsibility/Code/';
baseDir = [mfilepath filesep]; % should be Code directory
studies = {'Study 1 (behaviour only)','Study 2 (fMRI)'};
recpar_2st = {};
recpar_2st{1} = load([baseDir 'bin/compModels/Recovered model parameters for Behav study 10iter.mat']);
recpar_2st{2} = load([baseDir 'bin/compModels/Recovered model parameters for fMRI study 10iter.mat']);
actpar_2st = {};
actpar_2st{1} = readtable([baseDir 'csv/Behav - ResponsibilityRedux - fittedParameters.csv']);
actpar_2st{2} = readtable([baseDir 'csv/fMRI - ResponsibilityRedux - fittedParameters.csv']);
Npar = size(actpar_2st{1},2);

% Now, display
f2 = figure('Name','Comp model parameter recovery','position',[193 418 1313 472]);
fancyFig = 0;
fontsize = 16;
for st = 1:length(studies)
  recpar = recpar_2st{st}.recov_params_all;
  params = table2array(actpar_2st{st});
  
  % reorder, I want it to be: CR, EV, gamma, sRPE, self_pRPE
  newOrder = [1 2 4 3 5];
  recpar = recpar(:,newOrder,:);
  params = params(:,newOrder);
  paramNames = actpar_2st{st}.Properties.VariableNames(newOrder);
  paramNames{end} = 'social pRPE';
  Nrep = size(recpar,1);
  colors(2,:) = [.7 .7 .7];
  colors(1,:) = [.2 .2 .2]; %colors(2,:) + ([1 1 1] - colors(2,:))/2; % lighten up the color for the datapoints and the confidence intervals

  
  if fancyFig
    figure('name',['Parameter recovery for ' studies{st}],'position',[100 135 1261 806]);
    for p = 1:Npar
      subplot(Npar,1,p);
      violinDots(squeeze(recpar(:,p,:))); hold on; plot(params(:,p),'.k','markersize',16);
      ylabel(paramNames{p},'interpreter','none')
    end
  end
  
  figure(f2)
  for p = 1:Npar
    subplot(2,5,p+(st-1)*5);
    % plot(table2array(params(:,p)),squeeze(mean(recov_params_all(:,p,:))),'.') % x axis: real parameters across subjects, y axis: recovered parameters
    plot(repmat(params(:,p),1,Nrep)',squeeze(recpar(:,p,:)),'.') % x axis: real parameters across subjects, y axis: recovered parameters
%    regress_display(makerow(repmat(params(:,p),1,Nrep)'),makerow(squeeze(recpar(:,p,:))),...
    regress_display(params(:,p)',squeeze(mean(recpar(:,p,:))),...
      'inputNames',{paramNames{p},[paramNames{p} ' recov']},...
      'colors',colors,...
      'axesHandle',gca); % x axis: real parameters across subjects, y axis: recovered parameters
    set(gca,'fontsize',fontsize)
%     xlabel(paramNames{p},'interpreter','none')
%     ylabel([paramNames{p} ' recov'],'interpreter','none')
%    axis tight
%    if p == 3; title(studies{st},'fontsize',16); end
  title('')
  end
end
h = suptitle('Top row: Study 1; bottom row: Study 2');
set(h,'fontsize',fontsize+2)
legend('off')
if printFigFlag
  owd=pwd;
  cd([baseDir '/Figures/'])
  printfig
  printSVG
  cd(owd)
end