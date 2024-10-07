function [rho,sigma,Delta_rho,CI95] = fitCARAmodelToRiskyChoices2conds(valRiskyHigh,valRiskyLow,valSafe,chooseRisky,cond,plotFlag,figNameStr)
% function [rho,sigma,CI95] = fitCARAmodelToRiskyChoices2conds(valRiskyHigh,valRiskyLow,valSafe,chooseRisky,cond,plotFlag,figNameStr)
%
% This function fits an economic model of decision-making under risk to
% choice data (lottery vs. safe option).
%
% The model fitted is the constant absolute risk aversion model, which uses
% an exponential utility function. Thus, gain, mixed and loss trials need
% to be fitted separately.
%
% The utility function uses 1 parameter, rho, which indicates the level of
% risk aversion:
% rho > 0: risk averse
% rho = 0: risk neutral
% rho < 0: risk seeking
%
% The expected utility function allows to calculate the EU of a lottery
% based on the values of the low and high rewards, their probabilities, and
% the risk aversion parameter.
%
% From Holger's 2016 working paper:
% "A decision maker whose preferences can be represented by the subjective
% value function V (L;?) chooses A from {A,B} if ?V (A,B;?) > 0 and B if
% ?V (A,B;?) < 0. In this type of analysis, the subjective value V (L;?) is
% frequently set equal to the expected utility of the respective lottery."
%
% "The issue by using as the V -difference the difference between the
% lotteries? certainty equivalents. Formally, ?V(A,B;?) ?
% ?CE(A,B;?)=CE(A;?)?CE(B;?), (2) where CE(L;?) ? u?1?U(L;?);?? is the
% certainty equivalent of lottery L. U(L;?) denotes lottery L?s utility. We
% assume that it is given by expected utility..."
%
% In sum: the decision maker chooses the option with the highest certainty
% equivalent. CE of safe option = its value; CE of the risky option is
% determined by the EU of that option and the risk aversion parameter.
% Deviations from this optimal behaviour is described by sigma, Fechner
% noise
%
% So, we try to explain the subject's choices based on the difference
% between the CEs for the safe vs. risky options, with the risk aversion .
% The fit is implemented using a probit regression model with Fechner noise
% sigma, in which the Link function is the CDF of a Gaussian distribution.
%
% The fit is implemented using nonlinear least-squares regression
% (nlinfit.m).
%
% Without inputs, this code will show behavioural data of 3 theoretical
% participants.
%
% Problems: I have too few datapoints to only fit gain, loss or mixed
% trials.
% Solutions:
% 1) find another function that fits all trial types (prospect theory?)
% 2) fit data over all subjects at same time, with a subject random term
%
% Johannes Schultz and Holger Gerhardt
%
% Initial implementation summer 2020
% This version May 2021

demo = 0;
if nargin == 0
  close all;
  demo = 1;
end

if ~exist('plotFlag','var'); plotFlag = 1; end

% Which risk aversion framework to use?
useCRRA = 0; % if 0, uses CARA !! NOTE: CRRA does not work yet, can't simulate risk seeker & avoider data using my rho values

% Define utility functions for the estimation of the preference parameters of interest

% c is a commodity (here, money), and rho the risk aversion parameter

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

% CRRA (constant relative risk aversion)
% ---- Utility function:
 u_CRRA   = @(c, rho)   (rho ~= 1) .*    (c.^(1-rho) - 1) ./ (1 - rho - (rho == 1) .* 0.000001)  + ...
                        (rho == 1) .* log(c);
% ---- Expected utility function
EU_CRRA   = @(c1, c2, p1, rho)  p1 .* u_CRRA(c1, rho) + (1 - p1) .* u_CRRA(c2, rho);
% ---- Certainty equivalent under CRRA
% For regression as employed in von Gaudecker et al. (AER, 2011)
CE_CRRA   = @(EU, rho)  (rho ~= 1) .* ((1 - rho) .* EU + 1) .^ (1 ./ (1-rho-(rho==1).*0.000001)) + ...
                        (rho == 1) .*           exp(EU);
CE_CRRA_AiO = @(c1, c2, rho) CE_CRRA( EU_CRRA(c1, c2, 0.5, rho), rho );


% ---- Probit regression with Fechner noise sigma:
% Link function is CDF of Gaussian distribution
% Let CE1 be the risky option and CE2 the safe option
% sigma is noise
Link_probit = @(CE1, CE2, sigma)  1/2 .* (1 + erf((CE1 - CE2) ./ (sqrt(2) * sigma)));
if useCRRA == 1
  u_FCT       = u_CRRA;
  EU_FCT      = EU_CRRA;
  CE_FCT      = CE_CRRA;
  CE_FCT_AiO  = CE_CRRA_AiO;
else
  u_FCT       = u_CARA;
  EU_FCT      = EU_CARA;
  CE_FCT      = CE_CARA;
  CE_FCT_AiO  = CE_CARA_AiO;
end
Link_probit_AiO = @(cRiskyHigh, cRiskyLow, cSafe, cond, rho, sigma, Delta_rho) ...
  Link_probit(CE_FCT_AiO(cRiskyHigh, cRiskyLow, rho + Delta_rho .* cond), cSafe, sigma);

% Now a version with reorganised input arguments for use as MODELFUN for nlinfit:
% MODELFUN is a function, specified using @, that accepts two arguments, a
% coefficient vector and the array X, and returns a vector of fitted Y
% values.  BETA0 is a vector containing initial values for the coefficients.
% So I will set beta = [rho, sigma, Delta_rho] and X = [cRiskyHigh, cRiskyLow, cSafe, cond]
% Delta_rho is the difference in the risk aversion parameter between the
% self-other and the self-only condition (if cond==1 for self-other,
% Delta_rho > 0 indicates increased risk aversion in the self-other condition) 
Link_probit_AiO_nlinfit = @(beta,X) ...
  Link_probit_AiO(X(:,1), X(:,2), X(:,3), X(:,4), beta(1), beta(2), beta(3));

% For debugging, 2 versions with only 1 free parameter:
Link_probit_AiO_nlinfit_sigmaFixed = @(beta,X) ...
  Link_probit_AiO(X(:,1), X(:,2), X(:,3), X(:,4), beta(1), 0.5,       0); % sigma is fixed at 0.5
Link_probit_AiO_nlinfit_rhoFixed = @(beta,X) ...
  Link_probit_AiO(X(:,1), X(:,2), X(:,3), X(:,4), 0,       beta(1),   0); % rho is fixed at 0

Link = Link_probit_AiO_nlinfit; % needs as inputs beta = [rho, sigma] and X = [cRiskyHigh, cRiskyLow, cSafe]

% -------------------------------------------------------------------------
% Demo or fit real input data
if demo % create template subjects performing a set of decisions
  
  % create 3 subjects, plotted in different colours:
  if useCRRA
    rho       = [-0.1            0              0.1]'; %
  else % for CARA
    rho       = [-0.2            0              0.2]'; % !! for rho = 0.5, I get crazy certainty equivalents
  end
  subjNames = {'Risk seeker'; 'Risk neutral'; 'Risk avoider'};
  cols      = {'r';           'k';            'g'};
  T = table(rho,subjNames,cols);
  
  % Example lottery
  riskyHigh   = 8;
  riskyLow    = 4;
  probab      = 0.5;
  
  % Different values of safe alternative
  vSafe = riskyLow:.01:riskyHigh;
  
  % Theoretical noise: the dispersion (flatness) of the link function
  sigma = 0.5;
  
  % "real" noise: random decisions - in economics, they call this
  % "trembling"; using this would imply using another function than the
  % probit to model choices, so we don't use this anymore.
  noiseLevl = 0; % [0-1]
  
  % Create trials:
  Nt = 20000;
  simDatRiskyHigh = riskyHigh * ones(Nt,1);
  simDatRiskyLow  = riskyLow * ones(Nt,1);
  simDatValSafe   = linspace(riskyLow,riskyHigh,Nt)';
  
  % Create simulated choices of 3 template subjects:
  for s = 1:length(subjNames)
    chooseRisky{s} = Link_probit_AiO(simDatRiskyHigh, simDatRiskyLow, simDatValSafe, zeros(size(simDatRiskyHigh)), T.rho(s), sigma, 0) > 0.5;  % No noise: step function
    chooseRisky{s} = Link_probit_AiO(simDatRiskyHigh, simDatRiskyLow, simDatValSafe + randn(Nt,1) * sigma, zeros(size(simDatRiskyHigh)), T.rho(s), sigma, 0)  > 0.5;
      %abs(chooseRisky{s} + round(randn(Nt,1)*noiseLevl)); chooseRisky{s}(chooseRisky{s}>1) = 1; % add some noise - random additional decisions
  end
  T.chooseRisky = chooseRisky';
  
  % Create a figure to show the demo data
  figure('name','Demo CARA','position',[77         112        1077         581])
  
  % 1) show utility functions for risk seekers and avoiders
  subplot(2,4,1)
  c = -10:0.1:10; % Makes no sense for negative values
  for n = 1:numel(T.rho)
    plot(c,u_FCT(c,T.rho(n)),T.cols{n}); hold on
  end
  grid on
  legend(T.subjNames,'location','NorthWest')
  xlabel('Value')
  ylabel('Utility')
  title('Expon. utility fct for CARA')

  % 2) show expected utility of lottery for risk seekers and avoiders
  subplot(2,4,2)
  rhos_to_test = T.rho(1):0.001:T.rho(end);
  plot(rhos_to_test,EU_FCT(riskyHigh,riskyLow,probab,rhos_to_test))
  grid on
  title('EU of lottery as function of risk aversion \rho')
  ylabel(sprintf('EU of the lottery [%d,%d]; p=%.2f',riskyLow,riskyHigh,probab))
  xlabel('Rho (risk seeker <-> risk avoider)')
  
  % 3) show certainty equivalents for a given expected utility for risk seekers and avoiders
  subplot(2,4,3)
  EVlotteries = 2:12;
  for n = 1:numel(T.rho)
    plot(EVlotteries, CE_CARA_AiO(EVlotteries * 1.5, EVlotteries * 0.5, T.rho(n)), T.cols{n}); hold on
  end
  grid on
  title('CE for lotteries')
  legend(T.subjNames,'location','NorthWest')
  ylabel('CE')
  xlabel('EV of lottery')
  
  % 4) show P(safe) in the choice between the lottery as function of value
  % of safe option (for safe option: EV=SV=EU=CE) for our 3 template
  % subjects
  subplot(2,4,4); hold on
  for n = 1:numel(T.rho)
    plot(vSafe, 1-Link_probit_AiO(riskyHigh, riskyLow, vSafe, zeros(size(riskyHigh)), T.rho(n), sigma, 0), T.cols{n}); hold on
  end
  grid on
  title(sprintf('P(safe), in choice: [%d,%d]; p=%.2f vs. safe',riskyLow,riskyHigh,probab))
  legend(T.subjNames,'location','NorthWest')
  text(riskyLow*1.1,0.5,sprintf('sigma = %.2f',sigma))
  ylabel('P(safe)')
  xlabel('CE(safe)')
  
  % 5) Fit choices
  beta0 = [0 0.5 0]; % rho, sigma
%   sigmas_to_test = [.1:.1:2];
  nl_options = statset('MaxIter', 10000);
  
  for s = 1:length(subjNames)
%     for i = 1:length(sigmas_to_test)
%       beta01 = [beta0(1) sigmas_to_test(i)];
    [beta,resid,J,COVB] = nlinfit(...
      [simDatRiskyHigh,simDatRiskyLow,simDatValSafe,[zeros(size(simDatRiskyHigh)/2); ones(size(simDatRiskyHigh)/2); ] ],...
      T.chooseRisky{s},...
      Link, beta0, nl_options...
      );
    SSresid = resid'*resid;
    rho_est{s} = beta(1);    sigma_est{s} = beta(2);
    fprintf('%s: est. rho = %.3f; est. sigma = %.3f; SSresid = %.2f\n',T.subjNames{s},rho_est{s},sigma_est{s},SSresid)
%     end

    % 6) Parameter recovery: Show the similarity between original and
    % recovered parameters, and between original and fitted decisions
    chooseRisky_fit{s} = Link_probit_AiO(simDatRiskyHigh, simDatRiskyLow, simDatValSafe, rho_est{s}, sigma);
    
    % plot choices as function of CE safe, because lottery is always the
    % same - this plot should show non-overlapping lines for the 3 subjects
    subplot(2,4,5); hold on
    x = simDatValSafe;
    plot(x,1-T.chooseRisky{s},['o' T.cols{s}]); hold on
    hh(s) = plot(x,1-chooseRisky_fit{s},T.cols{s},'linewidth',2);
    plot([1 1]*((riskyHigh+riskyLow)/2),[0 1],[min(x) max(x)],[.5 .5],'-k')
    grid on
    ylabel('P(safe)')
    xlabel('CE(safe)')
    if s==1; title(sprintf('P(safe), in choice: [%d,%d]; p=%.2f vs. safe',riskyLow,riskyHigh,probab)); end

    % Now, plot risky choices as function of CElottery-CEsafe, where all lines should cross (0,0)
    subplot(2,4,5+s);
    % Get CEs:
    CElottery = CE_FCT_AiO(riskyHigh,riskyLow,rho_est{s});
    CEsafe    = simDatValSafe;
    x = CElottery-CEsafe;
    plot(x,T.chooseRisky{s},'o'); hold on
    plot(x,chooseRisky_fit{s},T.cols{s},'linewidth',2);
    plot([0 0],[0 1],[min(x) max(x)],[.5 .5],'-k')
    grid on
    ylabel('P(lottery)')
    xlabel('CE(lottery) - CE(safe)')
    legend({'Orig data','Fitted data'},'location','NorthWest')
    title(sprintf('real [%.3f %.3f]\nest. [%.3f %.3f]',T.rho(s),sigma,rho_est{s},sigma_est{s}))
  end % each subject
	legend(hh,T.subjNames,'location','NorthWest')

  T.chooseRisky_fit = chooseRisky_fit';

else % --------------- fit "real" input data ---------------------

  beta0       = [0 0.5 0]; % rho, sigma, Delta_rho
  nl_options  = statset('MaxIter', 10000);
  
  [beta,resid,J,COVB] = nlinfit(...
    [valRiskyHigh, valRiskyLow, valSafe, cond],...
    chooseRisky,...
    Link, beta0, nl_options...
    );
  rho = beta(1); sigma = beta(2); Delta_rho = beta(3);
  SSresid = resid'*resid;
  CI95 = nlparci(beta,resid,'covar',COVB); % 95% confidence intervals for beta, 1st row for rho, 2nd row for sigma
  fprintf('est. rho = %.3f; est. sigma = %.3f; est. Delta_rho = %.3f; SSresid = %.2f\n',rho,sigma,Delta_rho,SSresid)
  
  if plotFlag
    figure
    if exist('figNameStr','var'); set(gcf,'name',figNameStr); end
    % Show the similarity between original and fitted data
%     chooseRisky_fit = Link_probit_AiO(valRiskyHigh, valRiskyLow, valSafe, cond, rho, sigma, Delta_rho);
    
    CElottery       = CE_FCT_AiO(valRiskyHigh,valRiskyLow,rho+Delta_rho*cond);
    CEsafe          = valSafe;
    x               = CElottery-CEsafe;
    x_cont          = min(CElottery-CEsafe):0.1:max(CElottery-CEsafe);
    chooseRisky_fit = Link_probit(x_cont, zeros(size(x_cont)), sigma);
%     [x_sort,order]  = sort(x); chooseRisky_fit(order);
    
    plot(x,chooseRisky,'o'); hold on
%     plot(x,chooseRisky_fit,'*'); hold on
    plot(x_cont,chooseRisky_fit,'-','linewidth',2);
    plot([0 0],[0 1],'-k',[min(x) max(x)],[.5 .5],'-k')
    grid on
    ylabel('P(lottery)')
    xlabel('CE(lottery) - CE(safe)')
    legend({'Orig data','Fitted data'},'location','NorthWest')
    title(sprintf('rho: %.3f; sigma: %.3f; Delta rho: %.3f; SSresid: %.2f',rho,sigma,Delta_rho,SSresid))
  end % plot
end % demo or input data
