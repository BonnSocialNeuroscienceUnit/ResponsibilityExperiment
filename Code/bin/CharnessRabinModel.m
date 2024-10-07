function [alpha, beta, r2] = CharnessRabinModel(Thap)
% function [alpha, beta, r2] = CharnessRabinModel(Thap)
%
% Fits the Charness & Rabin (2002) inequality model to the happiness data.
% Input is table, with happiness and outcome values; outputs are the
% best-fitting parameters (grid search).
%
% In this analysis, I consider that happiness reflects utility, and I
% search for the parameters that minimze the difference between variations
% in utility and happiness.
%
% alpha = disadvantageous inequality parameter
% beta = advantageous inequality parameter
%
% Equation, search space and strategy from Smith & Krajbich 2018

alphaGrid = linspace(0,.5,11); % orig: 0:0.01:0.5
betaGrid = linspace(0,.75,11); % orig: 0:0.015:0.75

printIndivData = 0;

% For each subject, select only social trials in which both subject and
% partner obtained a reward, and where there was a happiness rating
for sub = 1:max(Thap.subject)
  idx = find(~isnan(Thap.happiness) & Thap.cond<3 & Thap.subject==sub);

  tt = 0;
  for t = 1:length(idx) % for each trial
    tt = tt + 1;
    xs = Thap.rewardSubj(idx(t)); % reward obtained by subject
    xp = Thap.rewardPart(idx(t)); % reward obtained by partner
    H(tt) = Thap.happiness(idx(t));
    if      xs > xp; s = 0; r = 1;
    elseif  xs < xp; s = 1; r = 0;
    else s = 0; r = 0;
    end
    for a = 1:length(alphaGrid)
      alphaT = alphaGrid(a);
      for b = 1:length(betaGrid)
        betaT = betaGrid(b);
        U(a,b,tt) = (1 - betaT*r - alphaT*s) * xs + (betaT*r + alphaT*s) * xp; % this is utility
        U2Hdiff(a,b,tt) = abs( U(a,b,tt) - H(tt) );
      end
    end
  end

  % Find best alpha and beta: those with smallest squared difference to happiness
  U2HdiffMean = nanmean(sqrt(U2Hdiff.^2),3);

  [a,b] = find(U2HdiffMean==min(min(U2HdiffMean)));

  alpha(sub) = mean(alphaGrid(a)) % in case there is more than 1 max, takes the average of all the best values
  beta(sub) = mean(betaGrid(b)) % in case there is more than 1 max, takes the average of all the best values

  % Display error matrix
  if printIndivData
    figure; imagesc(U2HdiffMean); xlabel('beta'); ylabel('alpha');
    set(gca,'xtick',1:length(alphaGrid),'ytick',1:length(betaGrid),'xticklabel',alphaGrid,'yticklabel',betaGrid)
  end

  % Evaluate whether the calculated utility fits the observed happiness
  % ratings
  Ufit = []; tt = 0;
  for t = 1:length(idx) % for each trial
    tt = tt + 1;
    xs = Thap.rewardSubj(idx(t)); % reward obtained by subject
    xp = Thap.rewardPart(idx(t)); % reward obtained by partner
    if isnan(xs); xs = 0; end
    if isnan(xp); xp = 0; end
    if      xs > xp; s = 0; r = 1; % subject got more than partner
    elseif  xs < xp; s = 1; r = 0; % subject got less than partner
    else s = 0; r = 0;
    end
    Ufit(tt) = (1 - beta(sub)*r - alpha(sub)*s) * xs + (beta(sub)*r + alpha(sub)*s) * xp; % the fitted utility
  end

  r2(sub) = (correl(Thap.happiness(idx),Ufit')^2);
end

figure; hist([alpha;beta]')