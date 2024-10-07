function out = exponUtilityFunction(beta,x)
% function exponUtilityFunction(beta,x)
% From Smith & Krajbich JEP Gen 2018

lambda = beta(1);
try b = beta(2);
catch b = 0; end

% x is the difference in objective value between the 2 options to choose
% from
% Lambda is the steepness of the curve

% This function is actually the same as logistic_js

% Example 
% figure
% x=-2:.1:2 % the differences in value between the 2 options
% for lambda=logspace(2,-3,100) % 100 values between 100 and 0.001
% plot(x,(1 + exp(-ll*(x))).^-1);
% grid on; hold on;
% end

out = (1 + exp(-lambda*(x+b))).^-1;