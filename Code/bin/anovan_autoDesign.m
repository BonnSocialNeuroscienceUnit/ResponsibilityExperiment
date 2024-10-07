function [T,str] = anovan_autoDesign(data,varargin)

% makes ANOVA design from data in a N-dimensional matrix, and passes on to
% anovan. Allows NaN in input matrix: ignores them.
%
% Varargin can be:
%     Parameter    Value
%     'alpha'      A value between 0 and 1 requesting 100*(1-alpha)%
%                  confidence bounds (default 0.05 for 95% confidence)
%     'continuous' A vector of indices indicating which grouping variables
%                  should be treated as continuous predictors rather than
%                  as categorical predictors
%     'display'    Either 'on' (the default) to display an anova table
%                  or 'off' to omit the display
%     'nested'     A matrix M of 0's and 1's specifying the nesting
%                  relationships among the grouping variables.  M(i,j)
%                  is 1 if variable i is nested in variable j.
%     'random'     A vector of indices indicating which grouping variables
%                  are random effects (all are fixed by default)
%     'sstype'     The type of sum of squares 1, 2, 3, or 'h' (default=3)
%     'varnames'   Grouping variables names in a character matrix or
%                  a cell array of strings, one per grouping variable
%                  (default names are 'X1', 'X2', ...)
%
%     'model'      The model to use, specified as one of the following:
%
%        'linear' to use only main effects of all factors (default)
%        'interaction' for main effects plus two-factor interactions
%        'full' to include interactions of all levels
%        an integer representing the maximum interaction order, for example
%           3 means main effects plus two- and three-factor interactions
%        a matrix of term definitions as accepted by the X2FX function,
%           but all entries must be 0 or 1 (no higher powers)
%
% J Schultz 23 Jun 2016 (Brexit vote day!!)

sizes = size(data);
Nfact = length(sizes);
sizes = [sizes 1]; % just for implementation

factors(:,1) = repmat(1:sizes(1),1,prod(sizes(2:end)))';
if Nfact>1;
  factors(:,2) = repmat(kron(1:sizes(2),ones(1,sizes(1))),1,prod(sizes(3:end)))';
  model = 2;
else; model = 1;
end
if Nfact>2;
  factors(:,3) = repmat(kron(1:sizes(3),ones(1,prod(sizes(1:2)))),1,prod(sizes(4:end)))';
end
if Nfact>3;
  factors(:,4) = repmat(kron(1:sizes(4),ones(1,prod(sizes(1:3)))),1,prod(sizes(5:end)))';
end
if Nfact>4;
  factors(:,5) = repmat(kron(1:sizes(5),ones(1,prod(sizes(1:4)))),1,prod(sizes(6:end)))';
end

d = data(:);

% find and remove NaN
ok = find(~isnan(d));
d = d(ok);
factors = factors(ok,:);

% read in options and format into anovan argument:
optionStr = [];
for v = 1:length(varargin)
  if ischar(varargin{v})
    optionStr = [optionStr ',''' varargin{v} ''''];
  elseif isnumeric(varargin{v}) % it's a number
    optionStr = [optionStr ',[' num2str(varargin{v}) ']'];
  elseif iscell(varargin{v}) % it's a cell: read each cell separately!
    optionStr = [optionStr ',{'];
    for c = 1:length(varargin{v})
      if c == 1
        optionStr = [optionStr '''' varargin{v}{c} ''''];
      else
        optionStr = [optionStr ',''' varargin{v}{c} ''''];
      end
    end
    optionStr = [optionStr '}'];
  end
end

% debug check:
% disp(optionStr)
% figure
% imagesc(normalise([d factors]))

% run it
% eval(['[P,T] = anovan(d,factors' optionStr ');'])
eval(['[P,T,stats] = anovan(d,factors' optionStr ',''model'',' num2str(model) ');'])

% print out partial eta squared (try to, works if 1st factor is a random effect)
try
  str = {};
  for tt = 2:size(T,1)-2
    eta_p2 = T{tt,2}/(T{tt,2}+T{end-1,2});
    DFdenom = T{tt,11};
    if round(DFdenom)==DFdenom
      str{end+1} = sprintf('%s: F(%d,%d) = %.2f, p = %.4f, eta_p^2 = %.2f',T{tt,1},T{tt,3},DFdenom,T{tt,6},T{tt,7},eta_p2);
    else
      str{end+1} = sprintf('%s: F(%d,%.1f) = %.2f, p = %.4f, eta_p^2 = %.2f',T{tt,1},T{tt,3},DFdenom,T{tt,6},T{tt,7},eta_p2);
    end
  end
  disp(char(str))
end

% Test if residuals are normally distributed
figure('position',[360   277   561   175],'name','Residuals')
[H,p,jbstat] = jbtest(stats.resid);
if H; str = 'Residuals not normally distributed'; else; str = 'Residuals normally distributed'; end
tit = sprintf('Jarque-Bera stat = %.1f, p = %.3f',jbstat,p);
subplot(1,2,1); plot(stats.resid); axis tight; ylabel('Residuals'); title(str)
subplot(1,2,2); qqplot(stats.resid); xlabel('Normal quantiles'); ylabel('Residuals'); title(tit)

% keyboard