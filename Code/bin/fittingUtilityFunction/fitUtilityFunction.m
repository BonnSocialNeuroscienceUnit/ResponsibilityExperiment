function fitUtilityFunction(decTab)

% from: https://www.mathworks.com/matlabcentral/answers/154879-using-fminsearch-to-estimate-utility-function-with-mle?s_tid=srchtitle

% I want to estimate the alpha (a) of a simple utility function of the form
% 
% EU(P,V) = P*V^a
%
% where P is probability of outcome, and V is magnitude of outcome. The
% response variable (riskChoice) is 1 for risky choice, and 0 for non-risky
% choice. I want to use maximum likelihood with a logit link-function to
% estimate the binary choice model:
% 
% 1/1(e^(w*(EUrisk^a-EUcertain^a))
%
% where EUrisk(cert) is the expected value of the risky(certain) option at
% each timepoint t, and w is some weight. From the above, I am interested
% in estimating w and a, using MLE.

riskChoice = decTab.chooseRisky;
% riskChoice = decTab.CARAfittedChoices; % fitted choices should make the function easier to fit but no! Why?
euCert = decTab.EVsafe;
euRisk = decTab.EVrisky;

% Binary choice model:
F = @(p) (1./(1 + exp(p(2)*(euCert.^p(1) - euRisk.^p(1)))));

logF1 = @(p) log(F(p)+eps);
logF2 = @(p) log(1-F(p)+eps);

% negative log-likelihood:
NLL = @(p) -sum( riskChoice.*logF1(p) + (1-riskChoice).*logF2(p) ) ;

Niter = 100;
initParam = zeros(Niter,numel(euRisk));
paramFrom = -10:0.1:10;
allOptWs=zeros(Niter,2); % a and w
modelEvi=zeros(Niter,1);
for r = 1:Niter
  a=rand*3;           % the exponent for value
  w=(rand-0.5)*20;    % a weight parameter (?)
  initWs = [a w];
  optWs = fminsearch(NLL,initWs);
  allOptWs(r,:) = optWs;
  modelEvi(r) = NLL(optWs);
end
[minModelEvi,minidx]=min(modelEvi);
optWs=allOptWs(minidx,:);

disp('Utility function estimated: EU(P,V) = P*V^a')
disp('Binary choice model: 1/1(e^(w*(EUrisk^a-EUcertain^a))')

a=optWs(1)
w=optWs(2)