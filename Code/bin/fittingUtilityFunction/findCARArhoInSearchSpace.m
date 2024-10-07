function rhoBest = findCARArhoInSearchSpace(riskyHigh,riskyLow,EVsafe,ChooseRisky,rhoToTest,plotFlag)
%
% function rhoBest = findCARArhoInSearchSpace(riskyHigh,riskyLow,EVsafe,ChooseRisky,rhoToTest);
%
% This function finds the value of the risk aversion parameter rho that,
% under the CARA model, best fits the observed choices between a lottery
% and a safe option.
%
% The input names are self explanatory.
%
% Johannes Schultz May 2021

if ~exist('plotFlag','var'); plotFlag = 1; end

if nargin == 0 % demo
  riskyHigh = 8*ones(50,1);
  riskyLow = 4*ones(50,1);
  EVsafe = linspace(riskyLow(1),riskyHigh(1),length(riskyHigh))';
  ChooseRisky = [ones(30,1); zeros(20,1)]; % a risk seeker!
  rhoToTest = .1:-0.005:-.1; [~,order] = sort(rhoToTest.^2); rhoToTest = rhoToTest(order); % here, I arrange the values of rho to test such that the value smallest to 0 are the most likely to be picked
end

% CARA (constant absolute risk aversion):
% ---- Utility function:
%  u_CARA   = @(c, rho) (1 - exp(-rho .* c)) ./ rho; % this is (1 - e^(-rho*c)) / rho, see https://en.wikipedia.org/wiki/Exponential_utility
u_CARA   = @(c, rho)   (rho ~= 0) .* (1 - exp(-rho .* c)) ./ (   rho + (rho == 0) .* 0.000001)  + ...
  (rho == 0) .*                  c;
% NOTE: rho is a constant that represents the degree of risk preference
% (rho>0 for risk aversion, rho=0 for risk-neutrality, and rho<0 for
% risk-seeking). c is a commodity, i.e. something the decision-maker likes

% ---- Expected utility function, which is the utility times the probability
% c1 and c2 are the high and low amounts of the risky option
EU_CARA   = @(c1, c2, p1, rho)  p1 .* u_CARA(c1, rho) + (1 - p1) .* u_CARA(c2, rho);

% ---- Certainty equivalent under CARA
% For regression as suggested by von Gaudecker et al. (AER, 2011, DOI: 10.1257/aer.101.2.664)
CE_CARA   = @(EU, rho)  (rho ~= 0) .*  log(1 - rho .* EU) ./ (   -rho + (rho == 0) .* 0.000001)  + ...
                        (rho == 0) .*                 EU;

% In this experiment, p1 = 0.5 for all trials, hence we do not need to
% carry this argument
CE_CARA_AiO = @(c1, c2, rho) CE_CARA( EU_CARA(c1, c2, 0.5, rho), rho );

% ----- ANALYSIS ----------------------------------------------------------
% now, create a set of certainty equivalents for the lotteries, with one
% set per rho value to test:
for i = 1:length(rhoToTest)
  CE_risky_CARA(:,i) = CE_CARA_AiO(riskyHigh,riskyLow,rhoToTest(i));
end

% then, compare the participant's choices with that of a decision-maker
% that chooses the option (lottery or safe) with the highest CE:
theoreticalChoices = CE_risky_CARA > EVsafe; % This fits very badly: theoreticalChoices = CE_risky_CARA - EVsafe;
% Calculate for each value of rho tested the deviation between the
% theoretical decision-maker and the actual choice:
diffBetweenChoices = mean((theoreticalChoices - ChooseRisky).^2);

% The best rho is the one for which the difference between theoretical
% choices and actual choices is the smallest:
[~,best] = min(diffBetweenChoices);
% rhoBest = rhoToTest(best);
rhoBest_idx = diffBetweenChoices==min(diffBetweenChoices);
rhoBest = mean(rhoToTest(rhoBest_idx));
if plotFlag
  figure;
  plot(rhoToTest,diffBetweenChoices,'.'); hold on
  plot(rhoToTest(rhoBest_idx),diffBetweenChoices(rhoBest_idx),'or','linewidth',2)
  xlabel('rho values to test')
  ylabel('Squared deviations between theoretical choices and real choices')
end
