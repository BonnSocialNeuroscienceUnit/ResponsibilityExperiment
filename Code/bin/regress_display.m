function [resultsall,mySTATSall,equationsall,betasall,hhall,y_covars] = regress_display(varargin)
% function [resultsall,mySTATSall,equationsall,betasall,hhall] = regress_display(x,y[,varargin])
% Uses regress.m to calculate simple regression between 2 input variables.
%
% Outputs a results string, a vector (or matrix) with [R2 F pval], the
% regression equations, the beta of the regression line, handles, and a
% vector representing the variance due to covariates.
%
% Can do robust regression by removing outliers (regressType = 'robust').
%
% Can regress out covariates specified as 'covariates'.
%
% Displays scatterplot with best-fitting regression line.
% Displays 95% confidence intervals if Pearson correlation is chosen.
% Can do different correlation types. Default is Pearson correlation.
% Outputs result string (R^2, F-value, p-value) also printed on the figure.
% Input variable names are optional (default is 'x' and 'y').
% If plotType is set to 2, displays 2 subplots, with in addition to the
% already-mentioned scatterplot a bicolor element-by-element display of the
% variables. If 0, does not plot.
% Colors is a 2X3 matrix defining the color of the data (1st row) and the
% fit line. Defaults are [0 0 1;0 0 1]: data and line blue, confidence area
% in lighter blue.

x = varargin{1};
y = varargin{2};

regressLineConfidenceAsArea = 1; % if 0, draw lines, if 1, paints area

% set defaults
correlType = 'Pearson'; % other possibilities: 'Spearman','Kendall','all',
regressType = 'OLS'; % other possibility: 'robust' to do robust regression; 'robust_nooutliers' to do robust regression but not show the outliers
showConfidence = 1; % other possibility: 0, then does not show confidence lines
R2Type = 'R^2'; % other possibility: 'R^2adj' for adjusted R-squared.
% inputNames = {'X','Y'};
plotType = 1; % others: 2 (see help), 0 = no plot.
% colors = [0 0 1;1 0 0]; % data blue, line red.
colors(2,:) = [1 .6 .03]; % nice warm orange
colors(1,:) = colors(2,:)+([1 1 1] - colors(2,:))/2; % lighten up the color for the datapoints and the confidence intervals
linewidth = 2;
markersize = 6;
threshCooksD = .5; % threshold to identify outliers using Cook's D
y_covars = []; % gets filled only if user gives a covariate

% set / get additional arguments, starting from 3rd input argument
% startOptArg = nargin+1;
for n = 1:nargin
  if ischar(varargin{n})
    % if startOptArg == nargin+1; startOptArg = n; end
    if strcmp(varargin{n},'colors')
      colors = varargin{n+1}; if size(colors,1) == 1; colors = repmat(colors,2,1);
      colors(1,:) = colors(2,:)+([1 1 1] - colors(2,:))/2; % lighten up the color for the datapoints and the confidence intervals
      end
    elseif strcmp(varargin{n},'correlType')
      correlType = varargin{n+1};
    elseif strcmp(varargin{n},'regressType')
      regressType = varargin{n+1};
    elseif strcmp(varargin{n},'showConfidence')
      showConfidence = varargin{n+1};
    elseif strcmp(varargin{n},'R^2Type')
      R2Type = varargin{n+1};
    elseif strcmp(varargin{n},'plotType')
      plotType = varargin{n+1};
    elseif strcmp(varargin{n},'inputNames')
      inputNames = varargin{n+1};
    elseif strcmp(varargin{n},'axesHandle')
      axesHandle = varargin{n+1};
    elseif strcmp(varargin{n},'linewidth')
      linewidth = varargin{n+1};
    elseif strcmp(varargin{n},'markersize')
      markersize = varargin{n+1};
    elseif strcmp(varargin{n},'covariates')
      covars = varargin{n+1};
    elseif strcmp(varargin{n},'titlePrefix')
      titlePrefix = varargin{n+1};
    end
  end
end

if exist('axesHandle','var'); plotType = 1; end % override plotType for plotting into predefined axes

% assess input, only take non-missing value pairs
if isempty(find(size(x)==1))
  error('Works only with 1-dimensional x axis variable')
end

if isempty(find(size(y)==1))
  disp('Columns are treated as separate variables')
else
  if size(y,2) ~= 1
    y = y';
  end
end

yall = y; x_orig = x;
try inputnamesall = inputNames; end
for v = 1:size(yall,2)
  y = yall(:,v); x = x_orig;
  try inputNames = inputnamesall([1 v+1]); end
  
  if size(x,2) ~= 1
    x = x';
  end
  
  % remove NaNs, and covariates if desired:
  idx = find(~isnan(x+y));
  if isempty('idx')
    error('No valid data: no pairs of values without missing data')
  end
  
  if exist('covars','var')
    y_in = y;
    if size(covars,2) == 1; temp = covars; else; temp = sum(covars')'; end
    idx = intersect(find(~isnan(x+y)),find(~isnan(temp)));
    y = y(idx,:);
    x = x(idx,:);
    covars = covars(idx,:);
    covars = [covars ones(size(covars,1),1)]; % add constant term
    b = regress(y,covars); % contributions of covariates to dependent variable
    %   [Q,R] = qr(covars,0);
    %   B = R\(Q'*y);
    y_covars = covars*b; % variance due to covariates
    y = y - y_covars; % subtract that from the dependent variables
    temp = NaN(size(y_in)); temp(idx) = y_covars; y_covars = temp; % y_covars in output should be the same length as the input
  else % just remove NaNs
    x = x(idx);
    y = y(idx);
  end
  
  % -------------------------------------------------------------------------
  % do regression, and stats
  if strfind(regressType,'robust') % use Cook's Distance to identify influential values & repeat regression without them
    nobs = length(find(~isnan(y)));
    X = [x ones(size(x,1),1)]; % add constant term
    [Q,R] = qr(X,0);
    B = R\(Q'*y);
    yhat = X*B;
    residuals = y - yhat;
    sse = norm(residuals)^2;    % sum of squared errors
    p = length(B);
    dfe = nobs-p;
    h = sum(abs(Q).^2,2);
    p = length(B);
    mse = sse./dfe;
    cooksd = abs(residuals).^2 .* (h./(1-h).^2)./(p*mse); % Cook's Distance
    try
      cooksd_interpreted = spm_Fcdf(cooksd,p-1,nobs-p); % Cook's D is looked up in F disribution, if it's <10-20 then it's not influential, but near 50 or more it is.
    catch
      cooksd_interpreted = fcdf(cooksd,p-1,nobs-p); % Cook's D is looked up in F disribution, if it's <10-20 then it's not influential, but near 50 or more it is.
    end
    % identify outliers:
    outliers = find(cooksd_interpreted>threshCooksD);
  else
    outliers = [];
  end
  
  if ~isempty(outliers) % then remove them
    regressTypeMessage = sprintf('-%d outliers',length(outliers));
    fprintf('Found %d outliers using Cook''s D, removing them and repeating regression\n',length(outliers))
    ok = find(cooksd_interpreted<=threshCooksD);
    x_out = x(outliers);
    y_out = y(outliers);
    x = x(ok);
    y = y(ok);
    X = [x ones(size(x,1),1)];
  end
  
  % calculate loads of things, taken from regstats.m:
  X = [x ones(size(x,1),1)]; % add constant term
  [Q,R] = qr(X,0);
  B = R\(Q'*y);
  yhat = X*B;
  
  residuals = y - yhat;
  nobs = length(find(~isnan(y)));
  p = length(B);
  
  dfe = nobs-p;
  dft = nobs-1;
  ybar = mean(y);
  
  sse = norm(residuals)^2;    % sum of squared errors
  ssr = norm(yhat - ybar)^2;  % regression sum of squares
  sst = norm(y - ybar)^2;     % total sum of squares;
  
  mse = sse./dfe;
  
  h = sum(abs(Q).^2,2);
  
  R2 = 1 - sse ./ sst;
  R2adj = 1 - (sse./sst)*(dft./dfe); % Adjusted R-square.
  
  F = (ssr/(p-1))/(sse/dfe);
  try
    pval = 1 - spm_Fcdf(F,p-1,nobs-p);   % Significance probability for regression
  catch
    pval = 1 - fcdf(F,p-1,nobs-p);   % Significance probability for regression
  end
  regressTypeMessage = regressType;
  
  switch R2Type
    case 'R^2'
      mySTATS = [R2 F ceil(pval*10000)/10000];
    case 'R^2adj'
      mySTATS = [R2adj F ceil(pval*10000)/10000];
  end
  
  switch correlType
    case 'Pearson'
      result = sprintf(['%s=%.2f, F(%d,%d)=%.2f, p=%.4f (%s)'],R2Type,mySTATS(1),p-1,nobs-p,mySTATS(2:3),regressTypeMessage);
    case 'Spearman'
      [rho,pval]=corr(x,y,'type',correlType);
      result = sprintf(['Rho=%.2f, p=%.4f'],rho,pval);
    case 'Kendall'
      [tau,pval]=corr(x,y,'type',correlType);
      result = sprintf(['Tau=%.2f, p=%.4f'],tau,pval);
    case 'all' % do all 3!
      result{1} = sprintf(['Regression: %s=%.2f, F=%.2f, p=%.4f (%s)'],R2Type,mySTATS,regressTypeMessage);
      [rho,pval]=corr(x,y,'type','Spearman');
      result{2} = sprintf(['Spearman''s Rho=%.2f, p=%.4f'],rho,pval);
      [tau,pval]=corr(x,y,'type','Kendall');
      result{3} = sprintf(['Kendall''s Tau=%.2f, p=%.4f'],tau,pval);
      result = char(result{:});
  end
  % result = sprintf(['R^2 = %.2f, F = %.2f, p < %.4f, RMSE = %.4f'],mySTATS,RMSE);
  
  % -------------------------------------------------------------------------
  % Stats values to calculate confidence interval for regression line
  % from here: http://shortrecipes.blogspot.de/2009/01/matlab-simple-linear-regression.html
  [xorder,order] = sort(x);
  n = length(x);
  alpha = 0.05; % if 0.05 -> 95% confidence intervals; this can be changed
  Syy=sum((y-mean(y)).^2);
  Sxy=sum((x-mean(x)).*(y-mean(y)));
  SSE=Syy-B(1)*Sxy;
  Syx=sqrt(SSE/(n-2)); % standard error of the estimate
  SSX=sum((x-mean(x)).^2);
  
  h = 1/n + ((xorder'-mean(xorder')).^2)/SSX;
  tn2 = tinv(1-alpha/2,n-2); %tn-2 - t critical
  interval = tn2.*Syx.*sqrt(h);
  
  % -------------------------------------------------------------------------
  if plotType ~= 0
    % Prepare figure, then plot
    if ~exist('axesHandle','var')
      hh(1) = figure('color',[1 1 1],'position',[576   580   441   420]);
      if size(yall,2) > 1
        subplot_optimal(size(yall,2),v)
        axesHandle = gca;
        hh(2) = axesHandle;
      else
        hh(2) = axes;
      end
    else
      if size(yall,2) > 1
        subplot_optimal(size(yall,2),v)
      end
    end
    
    if plotType == 2
      subplot(2,1,1); hold on
    else; hold on;
    end
    
    if exist('titlePrefix','var')
      result = [titlePrefix ': ' result];
    end
    title(result)
    % plot data
    pp(1) = plot(x,y,'o','color',colors(1,:),'markerfacecolor',colors(1,:),'markersize',markersize);
    if strcmpi(regressType,'robust') && length(outliers)>1
      pp(2) = plot(x_out,y_out,'x','color',colors(1,:),'markersize',markersize);
      legend(pp,'data','outlier')
    end
    % plot confidence on regression lines
    if (strcmpi(correlType,'Pearson') || strcmpi(correlType,'all')) && showConfidence
      if regressLineConfidenceAsArea
        hold on; hh(4:5) = errorarea(xorder,yhat(order),interval',colors(1,:),colors(1,:));
      else
        hh(4) = plot(xorder,yhat(order)+interval','color',colors(2,:),'linestyle','--');
        hh(5) = plot(xorder,yhat(order)-interval','color',colors(2,:),'linestyle','--');
      end
    end
    % plot regression line
    hh(3) = plot(xorder,yhat(order),'color',colors(2,:),'linewidth',linewidth);
    hh = [hh pp];
    
    % set(gca,'ylim',[(min(y)-.2) (max(y)+.2)],'xlim',[(min(x)-.2) (max(x)+.2)],'box','on');
    if exist('inputNames','var')
      xlabel(inputNames(1))%,'interpreter','none')
      ylabel(inputNames(2))%,'interpreter','none')
    else
      xlabel('X')
      ylabel('Y')
    end
    axis square, axis tight
    
    if plotType == 2
      subplot(2,1,2); hold on
      plot(y,'.r')
      plot(yhat)
      legend([inputNames(2) 'fitted values'])
    end
  end
  box on

%   % add R2 and p value on figure
%   text_relLoc(.1,.85,rr(15:22))
  
  % -------------------------------------------------------------------------
  equation = sprintf('yhat = %.3fx + %.3f',B(1),B(2));
  betas = B';
  
  %________________________________________________________________________
  % collect results for all columns of y
  resultsall{v} = result;
  mySTATSall{v} = mySTATS;
  equationsall{v} = equation;
  betasall{v} = betas;
  try hhall{v} = hh; end
end

if v == 1
  resultsall = resultsall{v};
  mySTATSall = mySTATSall{v};
  equationsall = equationsall{v};
  betasall = betasall{v};
  try hhall = hhall{v}; end
end

%plot(xorder,yhat(order))
%plot(y(order),x,'.r')
% keyboard
