function [rho,sigma,CI95] = fitCARAmodelToRiskyChoices(valRiskyHigh,valRiskyLow,valSafe,chooseRisky)

demo = 0;
if nargin == 0
  close all;
  demo = 1;
    
  % create 3 subjects, plotted in differnt colours:
  rho     = [-0.1            0              0.1]'; % !! for rho = 0.5, I get crazy certainty equivalents
  names   = {'Risk seeker'; 'Risk neutral'; 'Risk avoider'};
  cols    = {'r';           'k';            'g'};
  T = table(rho,names,cols);
  
  % Example lottery
  riskyHigh   = 6;
  riskyLow    = 3;
  probab      = 0.5;
  
  % Different values of safe alternative
  vSafe = riskyLow:.01:riskyHigh;
  
  % Noise; this is the inverse temperature of the curve
  sigma = 0.5;

  % Create simulated experimental data:
  Nt = 40;
  simDatRiskyHigh = riskyHigh * ones(Nt,1);
  simDatRiskyLow  = riskyLow * ones(Nt,1);
  simDatValSafe   = linspace(riskyLow,riskyHigh,Nt)';
  chooseRisky{1}  = zeros(Nt,1); chooseRisky{1}(1:round(Nt/2))   = 1; disp('real_rho = 0') % a risk-neutral subject
  chooseRisky{2}  = zeros(Nt,1); chooseRisky{2}(1:round(Nt/2.5)) = 1; disp('real_rho > 0')  % a risk-averse subject
  chooseRisky{3}  = zeros(Nt,1); chooseRisky{3}(1:round(Nt/1.5)) = 1; disp('real_rho < 0')   % a risk-seeking subject
%   chooseRisky   = [1 1 1 1 1 0 0 0 0 0]';   % a risk-neutral person
%   chooseRisky   = [1 1 1 1 1 1 1 0 0 0]'; % a risk-seeking person
%   chooseRisky   = [1 1 1 0 0 0 0 0 0 0]'; % a risk-avoiding person
  
  %   simDatRiskyHigh  = [30:10:130]';
  %   simDatRiskyLow   = [0:10:100]';
  %   simDatValSafe       = 65*ones(11,1);
  %   chooseRisky   = [0 0 0 1 0 0 0 0 0 1 1]';
  
end

% Define utility functions for the estimation of the preference parameters of interest

% c1 and c2 are the high and low amounts of the risky option

% CARA (constant absolute risk aversion):
% Utility function:
%  u_CARA   = @(c, rho) (1 - exp(-rho .* c)) ./ rho; % this is (1 - e^(-rho*c)) / rho, see https://en.wikipedia.org/wiki/Exponential_utility
u_CARA   = @(c, rho)   (rho ~= 0) .* (1 - exp(-rho .* c)) ./ (   rho + (rho == 0) .* 0.000001)  + ...
  (rho == 0) .*                  c;
% NOTE: rho is a constant that represents the degree of risk preference
% (rho>0 for risk aversion, rho=0 for risk-neutrality, and rho<0 for
% risk-seeking). c is a commodity, i.e. something the decision-maker likes

% Expected utility function, which is the utility times the probability
EU_CARA   = @(c1, c2, p1, rho)  p1 .* u_CARA(c1, rho) + (1 - p1) .* u_CARA(c2, rho);

% Certainty equivalent under CARA
% For regression as suggested by von Gaudecker et al. (AER, 2011, DOI: 10.1257/aer.101.2.664)
CE_CARA   = @(EU, rho)  (rho ~= 0) .*  log(1 - rho .* EU) ./ (   -rho + (rho == 0) .* 0.000001)  + ...
  (rho == 0) .*                 EU;
% In this experiment, p1 = 0.5 for all trials, hence we do not need to
% carry this argument
CE_CARA_AiO = @(c1, c2, rho) CE_CARA( EU_CARA(c1, c2, 0.5, rho), rho );

% Probit regression with Fechner noise sigma:
% Link function is CDF of Gaussian distribution
% Let CE1 be the risky option and CE2 the safe option
% sigma is noise
Link_probit = @(CE1, CE2, sigma)  1/2 .* (1 + erf((CE1 - CE2) ./ (sqrt(2) * sigma)));
Link_probit_AiO = @(cRiskyHigh, cRiskyLow, cSafe, rho, sigma) ...
  Link_probit(CE_CARA_AiO(cRiskyHigh, cRiskyLow, rho), cSafe, sigma);

% Now a version with reorganised input arguments for use as MODELFUN for nlinfit:
% MODELFUN is a function, specified using @, that accepts two arguments, a
% coefficient vector and the array X, and returns a vector of fitted Y
% values.  BETA0 is a vector containing initial values for the coefficients.
% So I will set beta = [rho, sigma] and X = [cRiskyHigh, cRiskyLow, cSafe]
Link_probit_AiO_nlinfit = @(beta,X) ...
  Link_probit_AiO(X(:,1), X(:,2), X(:,3), beta(1), beta(2));
Link_probit_AiO_nlinfit_sigmaFixed = @(beta,X) ...
  Link_probit_AiO(X(:,1), X(:,2), X(:,3), beta(1), 0.5); % sigma is fixed at 0.5
Link_probit_AiO_nlinfit_rhoFixed = @(beta,X) ...
  Link_probit_AiO(X(:,1), X(:,2), X(:,3), 0, beta(1)); % rho is fixed at 0

% Logit regression:
% Link function is CDF of logistic distribution
% Let CE1 be the risky option and CE2 the safe option
% sigma is noise
Link_logit  = @(CE1, CE2, sigma)   1 ./ (1 + exp(-(CE1 - CE2) ./ sigma));
Link_logit_AiO = @(cRiskyHigh, cRiskyLow, cSafe, rho, sigma) ...
  Link_logit( CE_CARA_AiO(cRiskyHigh, cRiskyLow, rho), cSafe, sigma );

Link = Link_probit_AiO; % needs as input cRiskyHigh, cRiskyLow, cSafe, rho, sigma

Link = Link_probit_AiO_nlinfit; % needs as inputs beta = [rho, sigma] and X = [cRiskyHigh, cRiskyLow, cSafe]

Link = Link_probit_AiO_nlinfit_sigmaFixed; % needs as inputs beta = [rho] and X = [cRiskyHigh, cRiskyLow, cSafe]

% TEST:
if demo
  figure('name','Demo utility function','position',[235   112   754   581])
  
  % 1) show utility functions for risk seekers and avoiders
  subplot(2,3,1)
  c = 0:10; % Makes no sense for negative values
  for n = 1:numel(T.rho)
    plot(c,u_CARA(c,T.rho(n)),T.cols{n}); hold on
  end
  grid on
  legend(T.names,'location','NorthWest')
  xlabel('Value')
  ylabel('Utility')
  title('Expon. utility fct for CARA')

  % 2) show expected utility of lottery for risk seekers and avoiders
  subplot(2,3,2)
  rhos_to_test = T.rho(1):0.001:T.rho(end);
  plot(rhos_to_test,EU_CARA(riskyHigh,riskyLow,probab,rhos_to_test))
  grid on
  title('EU of lottery as function of risk aversion')
  ylabel(sprintf('EU of the lottery [%d,%d]; p=%.2f',riskyLow,riskyHigh,probab))
  xlabel('Rho (risk seeker <-> risk avoider)')
  
  % 3) show certainty equivalents for a given expected utility for risk seekers and avoiders
  subplot(2,3,3)
  EVlotteries = 0:10;
  for n = 1:numel(T.rho)
    plot(EVlotteries,CE_CARA(EVlotteries,T.rho(n)),T.cols{n}); hold on
  end
  grid on
  title('CE for lotteries')
  legend(T.names,'location','NorthWest')
  ylabel('CE')
  xlabel('EV of lottery')
  
  % 4) show P(safe) in the choice between the lottery == [2,4] with p=0.5, as function of expected value of safe option for risk seekers and avoiders
  subplot(2,3,4); hold on
  for n = 1:numel(T.rho)
    plot(vSafe, 1-Link_probit_AiO(riskyHigh, riskyLow, vSafe, T.rho(n), sigma),T.cols{n}); hold on
  end
  grid on
  title(sprintf('P(safe), in choice: [%d,%d]; p=%.2f vs. safe',riskyLow,riskyHigh,probab))
  legend(T.names,'location','NorthWest')
  text(riskyLow*1.1,0.5,sprintf('sigma = %.2f',sigma))
  ylabel('P(safe)')
  xlabel('EV of safe option')
  
  % 5) Try to fit choices
  beta0 = [0 0.01]; % rho, sigma
  nl_options = statset('MaxIter', 10000);

%   disp(' ----- 1st fit: small N trials ------ ')
  Link = Link_probit_AiO_nlinfit_sigmaFixed; sigma = 0.5;
  beta0 = [0]; % beta = [rho sigma]
  [beta,resid,J,COVB] = nlinfit(...
    [simDatRiskyHigh,simDatRiskyLow,simDatValSafe],...
    chooseRisky,...
    Link, beta0, nl_options...
    );
  rho   = beta(1);
  fprintf('fitted rho = %.3f\n',rho)
%   sigma = beta(2)
  
  % 6) Parameter recovery: Show the similarity between original and recovered data
  chooseRisky_fit = Link_probit_AiO(simDatRiskyHigh, simDatRiskyLow, simDatValSafe, rho, sigma);
  
  subplot(2,3,5);
  plot(simDatValSafe,1-chooseRisky,'o'); hold on
  plot(simDatValSafe,1-chooseRisky_fit,'*r');
  ylabel('P(safe)')
  xlabel('EV of safe option')
  legend({'Orig data','Fitted data'},'location','NorthWest')
  title(sprintf('Fitted data, rho = %.3f',rho))

end

return
% keyboard

% % test: let's assume we have...
% rho = 0; % That's the risk aversion parameter: 0 for risk neutral, >0 for risk aversion, <0 for risk-seeking.
% sigma = 0.01; % That's the level of noise
% decisions_fitted = Link_probit_AiO( data_gain(:, 1), data_gain(:, 2), data_gain(:, 3), 0, 0.01 ); % This calls the Link function, which is currently Link_probit_AiO(CE_CARA_AiO(cRiskyHigh,cRiskyLow,rho),cSafe,sigma)
% % Recall: Arguments are (in this order) cRiskyHigh, cRiskyLow, cSafe rho, sigma
% 
% % decisions to fit = data(:,4); function to fit (Link) is
% % Link_probit_AiO, the independent variables are cRiskyHigh,
% % cRiskyLow and cSafe, and what we will estimate is rho (the risk
% % avoidance parameter) and sigma (the noise). These independent
% % variables and the default vaues of rho and sigma are inputs of
% % Link_probit_AiO function. To use nlinfit to estimate rho and sigma,
% % I created a new Link function that accepts rho and sigma as beta0,
% % and the independent variables as X.
% % nlinfit(X,y,MODELFUN,beta0,varargin)

% Define options and initial parameters
% nl_options = statset('MaxIter', 10000, 'Display', 'iter'); % shows each iteration
nl_options = statset('MaxIter', 10000);
% Recall: Link = Link_probit_AiO_nlinfit; % needs as inputs beta = [rho, sigma] and X = [cRiskyHigh, cRiskyLow, cSafe], to fit decisions stored for example in = data_gainOnly(:,4);

% Make 2nd set of simulation data, with larger N. Choice between a lottery
% in which the high and low rewards are 30 apart
Nt2 = 100; % larger N
valRiskyHigh  = linspace(30,130,Nt2)'; % the value of the risky
valRiskyLow   = linspace(0,100,Nt2)';
valSafe       = 65*ones(Nt2,1);
chooseRisky   = ones(Nt2,1); chooseRisky(1:round(Nt2/2))   = 0; % a risk-neutral subject
chooseRisky   = ones(Nt2,1); chooseRisky(1:round(Nt2/2.5))   = 0; % a risk-seeking subject
% chooseRisky   = ones(Nt2,1); chooseRisky(1:round(Nt2/1.5)) = 0; % a risk-averse subject
chooseRisky   = abs(chooseRisky + round(randn(Nt2,1)/3)); chooseRisky(chooseRisky>1) = 1; % add some noise - random additional decisions
  
% % =============== deactivated because did not work ========================
% disp('----- 2nd fit: large Ntrials ------')
% % Fit
% beta0 = [0 0.01]; % beta = [rho sigma]
% try [beta,resid,J,COVB] = nlinfit(...
%     [valRiskyHigh,valRiskyLow,valSafe],...
%     chooseRisky,...
%     Link, beta0, nl_options...
%     ); catch, beta = [NaN NaN]; end
% rho = beta(1)
% sigma = beta(2)
% CI95 = nlparci(beta,resid,'covar',COVB); % 95% confidence intervals for beta, 1st row for rho, 2nd row for sigma
% % =============== deactivated because did not work ========================

disp('----- 3rd fit: large N, only rho free to vary, sigma is fixed ------')
Link = Link_probit_AiO_nlinfit_sigmaFixed;
beta0 = [0]; % beta = [rho sigma]
try [beta,resid,J,COVB] = nlinfit(...
    [valRiskyHigh,valRiskyLow,valSafe],...
    chooseRisky,...
    Link, beta0, nl_options...
    ); catch, beta = [NaN NaN]; end
rho = beta(1)
% sigma = beta(2)

% =============== deactivated because did not work ========================
% disp('----- 4th fit: large N, rho is fixed, only sigma free to vary ------')
% Link = Link_probit_AiO_nlinfit_rhoFixed;
% beta0 = [0]; % beta = [rho sigma]
% try [beta,resid,J,COVB] = nlinfit(...
%     [valRiskyHigh,valRiskyLow,valSafe],...
%     chooseRisky,...
%     Link, beta0, nl_options...
%     ); catch, beta = [NaN NaN]; end
% % rho = beta(1)
% sigma = beta(1)
% =============== deactivated because did not work ========================




figure('name','CARA model','position',[223   341   703   314])

subplot(1,2,1)
EVrisky = (valRiskyHigh+valRiskyLow)/2;
% plot(EVrisky,EVrisky,'-*'); hold on
rhoVals = [-0.001 0 .001];
for v=1:length(rhoVals)
  for i=1:length(valRiskyHigh)
    EUrisky(i) = EU_CARA(valRiskyHigh(i),valRiskyLow(i),0.5,rhoVals(v));
  end
  hh(v) = plot(EVrisky,EUrisky); hold on
end
for i=1:length(valRiskyHigh)
  EUrisky(i) = EU_CARA(valRiskyHigh(i),valRiskyLow(i),0.5,rho);
end
hh(v+1) = plot(EVrisky,EUrisky,'k','linewidth',2); hold on

xlabel('Expected value of risky option')
ylabel('Expected utility of risky option')
legend(hh,num2str([rhoVals';rho]),'Location','NorthWest')
title('Utilities for different values of rho')

subplot(1,2,2)
plot(EUrisky-valSafe',chooseRisky','o')
xlabel('Expected utility of risky option - value of safe                                                                                                                                                                                                                                                                                                                                                                                                                                              ')
ylabel('P(choose risky option)')

% keyboard