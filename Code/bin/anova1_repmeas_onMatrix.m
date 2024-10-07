function [stats,message] = anova1_repmeas_onMatrix(data)
%
% function stats = anova1_repmeas(data)
%
% One-factor, within-subject repeated measures ANOVA.
% For designs with one within-subject factor.
%
% Parameters:
%    data       dependent variable (numeric), where:
%               ROWS are OBSERVATIONS
%               COLUMNS are FACTOR LEVELS or CONDITIONS
%
% Returns:
%    stats is a cell array with the usual ANOVA table:
%      Source / ss / df / ms / F / p
%
% Notes:
%    Program does not do any input validation, so it is up to you to make
%    sure that you have passed in the parameters in the correct form:
%
%       As data is a matrix, must have no missing values and 1 value per
%       condition/observation "cell".
%
% Johannes Schultz (2006.03.15), structure based on anova2_rp by Aaron Shurger

stats = cell(4,5);

F_levels = [1:size(data,2)];
Subjs = [1:size(data,1)];

a = length(F_levels); % # of levels in factor
n = length(Subjs); % # of subjects

INDS = cell(a,n); % this will hold arrays of indices
CELLS = cell(a,n); % this will hold the data for each subject X condition
MEANS = data'; % this will hold the means for each subj X condition


% compute subject, condition and grand means
subj_mean   = mean(MEANS,1);
cond_mean   = mean(MEANS,2)';
grand_mean  = mean(subj_mean);

% sums of squares
SStotal     = sum(sum(  (MEANS - repmat(grand_mean,a,n)).^2  ));
SScond      = n * sum(  (cond_mean - repmat(grand_mean,1,a)).^2  );
SSsubj      = a * sum(  (subj_mean - repmat(grand_mean,1,n)).^2  );
SSerror     = SStotal - SScond - SSsubj;

% degrees of freedom
DFcond      = a - 1;
DFsubj      = n - 1;
DFerror     = DFcond * DFsubj;

% mean of sums of squares
MScond      = SScond / DFcond;
MSsubj      = SSsubj / DFsubj;
MSerror     = SSerror / DFerror;

% F value
Fval        = MScond / MSerror;

% p value using fpdf
pval        = 1 - fcdf(Fval,DFcond,DFerror);

% partial eta squared:
partialEtaSq = SScond / (SScond + SSerror);

% return values
stats = {'Source','SS','df','MS','F','p','eta_p^2';...
    'Condition', SScond, DFcond, MScond, Fval, pval, partialEtaSq;...
    ['Subject'], SSsubj, DFsubj, MSsubj, [], [], [];...
    ['Error'], SSerror, DFerror, MSerror, [], [], [];...
    };
message = sprintf('F(%d,%d)=%.2f, p=%.4f, eta_p^2=%.3f',stats{2,3},stats{4,3},stats{2,5},stats{2,6},partialEtaSq);
if nargout == 0
    disp(message)
end

