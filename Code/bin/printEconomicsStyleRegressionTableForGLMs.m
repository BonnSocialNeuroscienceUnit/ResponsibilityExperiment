function printEconomicsStyleRegressionTableForGLMs(models, modelNames, titleStr, outputFilename)
%
% GLMs is a cell array of general linear models (from fitglm), or
% generalised linear mixed models (from fitglme).
% GLMnames is a cell array of the names of the models
% outputFilename is preferentially a .csv
%
% JS 9.12.2021

if ~exist('modelNames','var') || isempty(modelNames)
  for m = 1:length(models)
    modelNames{m} = ['Model ' num2str(m)];
  end
end

if ~exist('outputFilename','var')
  outputFilename = 'regressionOutput.csv';
end

precision = '2'; % N digits after comma for estimates
precision2 = '2'; % N digits after comma for standard errors of estimates

% First, gather all variable names, estimates, SEs and pvalues
% which will be the first table column
varNames = [];
estimates = [];
SEs = [];
pvals = [];
whichModel = [];
RsquaredO = [];
RsquaredA = [];
NumObs = [];
AIC = [];
BIC = [];
LogLik = [];

for m = 1:length(models)
  varNames  = [varNames; models{m}.CoefficientNames'];
  estimates = [estimates; models{m}.Coefficients.Estimate];
  SEs       = [SEs; models{m}.Coefficients.SE];
  pvals     = [pvals; models{m}.Coefficients.pValue];
  whichModel= [whichModel; m*ones(length(models{m}.CoefficientNames),1)];
  RsquaredO = [RsquaredO models{m}.Rsquared.Ordinary];
  RsquaredA = [RsquaredA models{m}.Rsquared.Adjusted];
  NumObs = [NumObs models{m}.NumObservations];
  AIC = [AIC models{m}.ModelCriterion.AIC];
  BIC = [BIC models{m}.ModelCriterion.BIC];
  LogLik = [LogLik models{m}.LogLikelihood];
end

% sort variable names to form first column and index to gather data for
% next columns
[C,~,IC] = unique(varNames,'stable');

% gather information in final table form
tabPval = NaN(length(C),m);
tabEst = NaN(length(C),m);
tabSE = NaN(length(C),m);
for v = 1:length(estimates)
  tabEst(IC(v),whichModel(v)) = estimates(v);
  tabSE(IC(v),whichModel(v)) = SEs(v);
  tabPval(IC(v),whichModel(v)) = pvals(v);
end

% create output cells: add stars to estimates, replace NaNs with ''
for i = 1:size(tabEst,1)
  for j = 1:size(tabEst,2)
    star = stars(tabPval(i,j));
    if ~isnan(tabEst(i,j))
      tabEstStr{i,j} = sprintf(['%.' precision 'f%s'],tabEst(i,j),star{1});
      tabSEStr{i,j} = sprintf(['(%.' precision2 'f)'],tabSE(i,j));
    else
      tabEstStr{i,j} = '';
      tabSEStr{i,j} = '';
    end
  end
end

% print to file
fid = fopen(outputFilename,'w');

if exist('titleStr','var') && ~isempty(titleStr)
  fprintf(fid,'%s\n',titleStr);
end

fprintf(fid,[' ' repmat(', %s',1,m) '\n'],modelNames{:});
for r = 1:size(tabEst,1)
  fprintf(fid,['%s ' repmat(', %s',1,m) '\n'],C{r}, tabEstStr{r,:});
  fprintf(fid,['%s ' repmat(', %s',1,m) '\n'],'           ', tabSEStr{r,:});
end
fprintf(fid,['R2 (ord) ' repmat(', %.3f',1,m) '\n'],RsquaredO);
fprintf(fid,['R2 (adj) ' repmat(', %.3f',1,m) '\n'],RsquaredA);
fprintf(fid,['AIC ' repmat(', %.0f',1,m) '\n'],AIC);
fprintf(fid,['BIC ' repmat(', %.0f',1,m) '\n'],BIC);
fprintf(fid,['LogLikelihood ' repmat(', %.0f',1,m) '\n'],LogLik);
fprintf(fid,['N ' repmat(', %.d',1,m) '\n'],NumObs);
fclose(fid);


return
fprintf(1,[' ' repmat('\t %s',1,m) '\n'],GLMEnames{:});
for r = 1:size(tabEst,1)
  fprintf(1,['%s ' repmat('\t %s',1,m) '\n'],C{r}, tabEstStr{r,:});
  fprintf(1,['%s ' repmat('\t %s',1,m) '\n'],'           ', tabSEStr{r,:});
end

% ---- second approach: make combined table

% add estimates and SEs to table
bigTab = [];
bigTab(1:2:length(IA)*2,:) = tabEst;
bigTab(2:2:length(IA)*2,:) = tabSE;

% adapt 1st column by adding empty cells
bigTabC(1:2:length(IA)*2,:) = C;
bigTabC(2:2:length(IA)*2,1) = {''};

% simple display
[char(bigTabC) num2str(bigTab)]

