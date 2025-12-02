function [result] = fit_happy_model_responsibility_redux(mtx,happyscore,constant)

%assumes that each row of the cell array is a different subject

rptfit = 5; %repeat N times jittering the betas
tempf=fieldnames(mtx);
for n=1:length(tempf), eval(sprintf('mtx.%s=double(mtx.%s);',tempf{n},tempf{n})); end;
data = mtx;
data.happyscore = double(happyscore);

%result = data;
if isfield(mtx,'inx')
    inx = mtx.inx;
else
    inx = [0.5 0.5 0.5 0.5    0]; % dummy parameter values, 1 per regressor
end
options = optimset('Display','off','MaxIter',100000,'TolFun',1e-10,'TolX',1e-10,...
    'DiffMaxChange',1e-2,'DiffMinChange',1e-4,'MaxFunEvals',10000,'LargeScale','off');
warning off; %display,iter to see outputs
lb = [-5 -5 -5 0.01   -1]; % lower bounds of parameter values, 1 per regressor
ub = [5 5 5 1          1]; % upper bounds of parameter values, 1 per regressor

if exist('constant','var') && ~isnan(constant),
    inx = [inx mean(happyscore)];
    lb = [lb 0];
    ub = [ub 1];
    if isfield(mtx,'const') %use previous fit constant
        inx(end) = mtx.const;
    end;
end;

result.modelName = 'ResponsibilityRedux';
result.inx = double(inx);
dof = length(inx);
%result.options = options;
result.lb  = lb;
result.ub  = ub;
b=inx;
%b = fmincon(@model_param, inx, [],[],[],[],lb,ub,[], options, data);
result.b = b;
[lse, happypred, r2, r2adj] = model_param(b, data);

for n=1:rptfit, %jitter the betas
    b = fmincon(@model_param, double(b+0.1*randn(size(b))), [],[],[],[],lb,ub,[], options, data);
    lse2 = model_param(b, data);
    if lse2<lse,
        result.b = b;
        [lse, happypred, r2, r2adj] = model_param(b, data);
    end;
end;

result.happypred = happypred;
result.r2 = r2;
result.r2adj = r2adj;

% calculate bic, aic and log likelihood
k = length(b); % 𝑘: The number of estimated parameters in the model.
n = length(happypred); % n: The number of data points or observations used in the model.
mse = lse/n; % mse: mean squared error

result.bic = n*log(mse) + k*log(n);
result.aic = n*log(mse) + 2*k;

% Now I want the log-likelihood of the model fit, so I can run a
% likelihood ratio test: LRT = -2 * ln(L0/L1) = -2 * (logL0 - logL1). The LRT
% value often follows a chi-square distribution. The degrees of freedom are
% the difference in the number of parameters between the two models.

%
% In the apparently standard formulae:
% bic = -2*log(L) + k*log(n);
% aic = -2*log(L) + 2*k;
%
% L is the maximized value of the likelihood function for the model. This
%  measures how well the model fits the data. In our case, we run a model
%  that assumes a normal (Gaussian) distribution with constant variance for
%  the errors. Under these assumtions, the negative log-likelihood function
%  is proportional to the mean squared error (MSE).
%
% Comparing the formulae for AIC and BIC, it follows that:
%  -2*log(L) = n*log(mse)
% Thus we should have
%  log(L) = n*log(mse) / -2
%
% The check worked: I got the same bic and aic using the two versions of
% these formulae.
%  check_bic = -2*logL + k*log(n);
%  check_aic = -2*logL + 2*k;
result.logL = n*log(mse) / -2;

result.paramNames = {'CR','EV','sRPE','gamma','social_pRPE'};
if exist('constant','var') && ~isnan(constant)
  result.paramNames = [result.paramNames, {'constant'}];
end


function [lse, happypred, happyr2, happyr2adj] =  model_param(x, data)

a = x(1);   %certain
b = x(2);   %ev
c = x(3);   %rpe
tau = x(4); %decay constant
try s1 = x(5);  %partner RPE resulting from self decisions, regressor called "self"
end
if length(x)==6,
    const = x(6);
else
    const = 0; %if z-scored
end;

%needs all blocks to be the same size
decayvec = tau.^[0:size(data.certainmtx,2)-1]; decayvec = decayvec(:);
theychose = data.theychosemtx; youchose = data.youchosemtx; cert = data.certainmtx; ev = data.evmtx; rpe = data.rpemtx; otherrpe = data.otherrpemtx; dec = decayvec;

otherrpeactive = youchose.*data.otherrpemtx; %partner RPE resulting from self decisions, regressor called "self"
otherrpepassive = theychose.*data.otherrpemtx; %partner RPE resulting from partner decisions, regressor called "partner"

% HERE THE MODEL IS DEFINED
if exist('s1','var')
  happypred = a*cert*dec + b*ev*dec + c*rpe*dec + ...
    s1*otherrpeactive*dec + const;
else
  happypred = a*cert*dec + b*ev*dec + c*rpe*dec + ...
    const;
end

try lse = sum((data.happyscore-happypred).^2); %sum least-squares error
catch lse = sum((data.happyscore-happypred(2:2:end)).^2); %sum least-squares error for fMRI regressors
end

meanhappy = mean(data.happyscore);
re=sum((data.happyscore-meanhappy).^2);
happyr2 = 1-lse/re;

nobs = length(find(~isnan(data.happyscore)));
p = 5;

happyr2adj = 1 - (lse/re) * ((nobs-1)/(nobs-p));
