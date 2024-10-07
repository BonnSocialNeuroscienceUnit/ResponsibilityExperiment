function [yfit,beta] = fitExponUtilityFunction(x,y,plotFlag)
% function [yfit,beta] = fitExponUtilityFunction(x,y,plotFlag)
%
% This function fits the logistic function implemented in the code
% exponUtilityFunction.m using ordinary least sqaures.

startBeta = [1 1];
func = 'exponUtilityFunction'; beta = lsqcurvefit(func,startBeta,x,y);
yfit = exponUtilityFunction(beta,x);
if plotFlag
  figure; plot(x,yfit,'o-');
end
