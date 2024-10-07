function [d,r] = cohens_d(data1,data2,pairedSamplesFlag)
% function [d,r] = cohens_d(data1[,data2,pairedSamplesFlag])
%
% Calculates a standardized measure of effect size for t-tests:
% Cohen's d and r for 1-sample and 2-sample t-tests.
% INPUTS:
%   data1:  mandatory, vector for 1st data sample.
%   data2:  if scalar, is used as reference nanmean for 1-sample test. If
%           vector, used as data for 2nd sample. If absent, 1-sample t
%           against 0 is done.
%   pairedSamplesFlag: only used for 2-sample tests with equal N: 1 for
%           paired samples t, 0 for independent samples. If absent, data2
%           exists and is same length as data1, it's requested from user.
%
% OUTPUTS:  
%   d:      Cohen's d. 0.2 is small effect, 0.5 medium effect, 0.8+ is
%           large effect.
%   r:      "coefficient of intraclass correlation" - the good old r value.
%
% johannes - march 2008
% Refs: 
% - Jacob Cohen, "Statistical power analysis for the behavioral
% sciences." Book, 1988.
% - Cohen, J. [...] A power primer. Psych Bull, 112:155-159 (1992).
% - www.en.Wikipedia.org/wiki/Effect_size
% Implementation: http://www2.chass.mcsu.edu/garson/PA765/anova.htm#cohensd

if ~exist('data2','var')
    str = '1-sample t-test, comparing nanmean against 0';
    d = abs( nanmean(data1) ) ./ nanstd(data1);
else
    if length(data2) == 1
        str = sprintf('1-sample t-test, comparing nanmean against %d',data2);
        d = abs( ( nanmean(data1) - data2 ) ) ./ nanstd(data1);
    else
        if length(data2) ~= length(data1)
            pairedSamplesFlag = 0;
        elseif length(data2) == length(data1)
            if ~exist('pairedSamplesFlag','var')
                pairedSamplesFlag = input('Paired t-test, or not? [1/0]: ');
            end
        end
        if pairedSamplesFlag == 1
            str = 'paired samples t-test';
            data1 = data2-data1;
            d = abs(nanmean(data1)) ./ nanstd(data1);
        elseif pairedSamplesFlag == 0
            str = 'independent samples t-test';
            d = abs( ( nanmean(data1) - nanmean(data2) ) ) ./ ...
                sqrt( ( (length(data1)-1)*nanstd(data1).^2 + (length(data2)-1)*nanstd(data2).^2 ) ...
                / (length(data1)+length(data2)-2) );
        end
    end
end

r = d ./ (d.^2+4).^.5;

if ~nargout
    fprintf('%s: d=%.4f, r=%.4f\n',str,d,r)
end
