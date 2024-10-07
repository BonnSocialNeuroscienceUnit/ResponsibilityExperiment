function CI = ci_t(data,percentile)
% function CI = ci_t(data[,percentile])
%
% Returns confidence intervals for the t-distribution for normally distributed data.
%
% Formula: CI = m +/- t_crit * (std / sqrt(N))
% where m is the mean of the data sample, t_crit is the critical value of
% the t distribution, std is the standard deviation of the sample, and N is
% the sample size.
%
% if data is a MxNxS matrix, operates along columns, and gives out a 2xNxS
% matrix of confidence intervals.

if ~exist('percentile','var')
    percentile = 95;
end

t_crit_location = [(100-percentile)/2 percentile + (100-percentile)/2]/100; % [0.025 0.975] for 95% confidence intervals

if any(size(data)==1)
    data = data(:);
end
for c = 1:size(data,2) % each column
  for s = 1:size(data,3) % each slice of a 3D matrix
    dat     = data(:,c,s);
    dat     = dat(~isnan(dat));
    SEM     = nanstd(dat)/sqrt(length(dat));               % Standard Error
    ts     = tinv(t_crit_location, length(dat)-1);      % T-Score
    CI(:,c,s) = nanmean(dat) + ts*SEM;                      % Confidence Intervals
  end
end
