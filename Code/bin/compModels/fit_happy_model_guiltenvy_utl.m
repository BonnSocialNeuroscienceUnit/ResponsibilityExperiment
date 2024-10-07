function [result] = fit_happy_model_nogain_guiltenvy_utl(mtx,happyscore,constant)

%assumes that each row of the cell array is a different subject

rptfit = 100; %repeat N times jittering the betas
tempf=fieldnames(mtx);
for n=1:length(tempf), eval(sprintf('mtx.%s=double(mtx.%s);',tempf{n},tempf{n})); end;
data = mtx;
data.happyscore = double(happyscore);

%result = data;
if isfield(mtx,'inx'),
    inx = mtx.inx;
else
    inx = [0.5 0.5 0.5 0.5   0 0 0.5]; %dummies 
end
options = optimset('Display','off','MaxIter',100000,'TolFun',1e-10,'TolX',1e-10,...
    'DiffMaxChange',1e-2,'DiffMinChange',1e-4,'MaxFunEvals',10000,'LargeScale','off');
warning off; %display,iter to see outputs
lb = [-5 -5 -5 0.01   -1 -1 0.005];
ub = [5 5 5 1          1  1 100];

if exist('constant','var') && ~isnan(constant),
    inx = [inx mean(happyscore)];
    lb = [lb 0];
    ub = [ub 1];
    if isfield(mtx,'const'), %use previous fit constant
        inx(end) = mtx.const;
    end;
end;

result.modelName = 'GuiltEnvy';
result.inx = double(inx);
dof = length(inx);
%result.options = options;
result.lb  = lb;
result.ub  = ub;
b=inx;
%b = fmincon(@model_param, inx, [],[],[],[],lb,ub,[], options, data);
result.b = b;
[lse, happypred, r2] = model_param(b, data);

for n=1:rptfit, %jitter the betas
    b = fmincon(@model_param, double(b+0.1*randn(size(b))), [],[],[],[],lb,ub,[], options, data);
    lse2 = model_param(b, data);
    if lse2<lse,
        result.b = b;
        [lse, happypred, r2] = model_param(b, data);
    end;
end;

result.happypred = happypred;
result.r2 = r2;
result.bic = length(happypred)*log(lse/length(happypred)) + length(b)*log(length(happypred));
result.aic = length(happypred)*log(lse/length(happypred)) + 2*length(b);
result.paramNames = {'certain','ev','rpe','decay','guilt','envy','Utl'};
if exist('constant','var') && ~isnan(constant)
  result.paramNames = [result.paramNames, {'constant'}];
end


function [lse, happypred, happyr2] =  model_param(x, data)

a = x(1);   %certain
b = x(2);   %ev
c = x(3);   %rpe
tau = x(4); %decay constant
try s1 = x(5);  %guilt
  s2 = x(6);  %envy
end
utl = x(7); %decay constant
if length(x)==8,
    const = x(8);
else
    const = 0; %if z-scored
end;

%needs all blocks to be the same size
decayvec = tau.^[0:size(data.certainmtx,2)-1]; decayvec = decayvec(:);
happypred = data.happyscore;
theychose = data.theychosemtx; youchose = data.youchosemtx; cert = data.certainmtx; ev = data.evmtx; rpe = data.rpemtx; dec = decayvec;
%ywtw = data.notguiltmtx;
ywtl = data.guiltmtx; %double(data.guiltmtx>0); %or leave as parametric - how much they got less than you
yltw = data.envymtx; %double(data.envymtx>0); %or leave as parametric - how much they got more than you
%yltl = data.notenvymtx;
utl2=(utl.^2);
range=14;

term1x=cert*dec;
term1x2=zeros(size(term1x,1),1);
for res_loi=1:size(term1x2,1)
if term1x(res_loi,1)<0
% term1x2(res_loi,1)=-log(-term1x(res_loi,1))- 2.*utl.*((range/2)+(exp(term1x(res_loi,1))/2));%Minus Utility function
% term1x2(res_loi,1)=real(term1x2(res_loi,1));
 term1x2(res_loi,1)=a.*((range./(1+exp(-(3/4).*(utl).*(term1x(res_loi,1)))))-(range/2));
if term1x2(res_loi,1)< -(range/2)
term1x2(res_loi,1)= -(range/2);
elseif term1x2(res_loi,1)==-Inf
  term1x2(res_loi,1)= -(range/2);
end
elseif term1x(res_loi,1)>0
% term1x2(res_loi,1)=log(term1x(res_loi,1))+ utl.*((range/2)+(exp(-term1x(res_loi,1))/2)); %Positive Utility function
 term1x2(res_loi,1)=a.*((range./(1+exp(-(3/12).*(utl).*(term1x(res_loi,1)))))-(range/2));
if term1x2(res_loi,1)> (range/2)
    term1x2(res_loi,1)= (range/2);
end
end
end
term2x=ev*dec;
term2x2=zeros(size(term2x,1),1);
for res_loi=1:size(term2x2,1)
if term2x(res_loi,1)<0
% term1x2(res_loi,1)=-log(-term1x(res_loi,1))- 2.*utl.*((range/2)+(exp(term1x(res_loi,1))/2));%Minus Utility function
% term1x2(res_loi,1)=real(term1x2(res_loi,1));
 term2x2(res_loi,1)=b.*((range./(1+exp(-(3/4).*(utl).*(term2x(res_loi,1)))))-(range/2));
if term2x2(res_loi,1)< -(range/2)
term2x2(res_loi,1)= -(range/2);
elseif term2x2(res_loi,1)==-Inf
  term2x2(res_loi,1)= -(range/2);
end
elseif term2x(res_loi,1)>0
% term1x2(res_loi,1)=log(term1x(res_loi,1))+ utl.*((range/2)+(exp(-term1x(res_loi,1))/2)); %Positive Utility function
 term2x2(res_loi,1)=b.*((range./(1+exp(-(3/12).*(utl).*(term2x(res_loi,1)))))-(range/2));
if term2x2(res_loi,1)> (range/2)
    term2x2(res_loi,1)= (range/2);
end
end
end
term3x=rpe*dec;
term3x2=zeros(size(term3x,1),1);
for res_loi=1:size(term3x2,1)
if term3x(res_loi,1)<0
% term1x2(res_loi,1)=-log(-term1x(res_loi,1))- 2.*utl.*((range/2)+(exp(term1x(res_loi,1))/2));%Minus Utility function
% term1x2(res_loi,1)=real(term1x2(res_loi,1));
 term3x2(res_loi,1)=c.*((range./(1+exp(-(3/4).*(utl).*(term3x(res_loi,1)))))-(range/2));
if term3x2(res_loi,1)< -(range/2)
term3x2(res_loi,1)= -(range/2);
elseif term3x2(res_loi,1)==-Inf
  term3x2(res_loi,1)= -(range/2);
end
elseif term3x(res_loi,1)>0
% term1x2(res_loi,1)=log(term1x(res_loi,1))+ utl.*((range/2)+(exp(-term1x(res_loi,1))/2)); %Positive Utility function
 term3x2(res_loi,1)=c.*((range./(1+exp(-(3/12).*(utl).*(term3x(res_loi,1)))))-(range/2));
if term3x2(res_loi,1)> (range/2)
    term3x2(res_loi,1)= (range/2);
end
end
end
term4x=ywtl*dec;
term4x2=zeros(size(term4x,1),1);
for res_loi=1:size(term4x2,1)
if term4x(res_loi,1)<0
% term1x2(res_loi,1)=-log(-term1x(res_loi,1))- 2.*utl.*((range/2)+(exp(term1x(res_loi,1))/2));%Minus Utility function
% term1x2(res_loi,1)=real(term1x2(res_loi,1));
 term4x2(res_loi,1)=s1.*((range./(1+exp(-(3/4).*(utl).*(term4x(res_loi,1)))))-(range/2));
if term4x2(res_loi,1)< -(range/2)
term4x2(res_loi,1)= -(range/2);
elseif term4x2(res_loi,1)==-Inf
  term4x2(res_loi,1)= -(range/2);
end
elseif term4x(res_loi,1)>0
% term1x2(res_loi,1)=log(term1x(res_loi,1))+ utl.*((range/2)+(exp(-term1x(res_loi,1))/2)); %Positive Utility function
 term4x2(res_loi,1)=s1.*((range./(1+exp(-(3/12).*(utl).*(term4x(res_loi,1)))))-(range/2));
if term4x2(res_loi,1)> (range/2)
    term4x2(res_loi,1)= (range/2);
end
end
end
    term5x=yltw*dec;
term5x2=zeros(size(term5x,1),1);
for res_loi=1:size(term5x2,1)
if term5x(res_loi,1)<0
% term1x2(res_loi,1)=-log(-term1x(res_loi,1))- 2.*utl.*((range/2)+(exp(term1x(res_loi,1))/2));%Minus Utility function
% term1x2(res_loi,1)=real(term1x2(res_loi,1));
 term5x2(res_loi,1)=s2.*((range./(1+exp(-(3/4).*(utl).*(term5x(res_loi,1)))))-(range/2));
if term5x2(res_loi,1)< -(range/2)
term5x2(res_loi,1)= -(range/2);
elseif term5x2(res_loi,1)==-Inf
  term5x2(res_loi,1)= -(range/2);
end
elseif term5x(res_loi,1)>0
% term1x2(res_loi,1)=log(term1x(res_loi,1))+ utl.*((range/2)+(exp(-term1x(res_loi,1))/2)); %Positive Utility function
 term5x2(res_loi,1)=s2.*((range./(1+exp(-(3/12).*(utl).*(term5x(res_loi,1)))))-(range/2));
if term5x2(res_loi,1)> (range/2)
    term5x2(res_loi,1)= (range/2);
end
end
end
if exist('s1','var')
  happypred = (term1x2 + term2x2 +term3x2 +term4x2 +term5x2)./5;
else
  happypred = (term1x2 + term2x2 +term3x2)./3;
end

lse = sum((data.happyscore-happypred).^2); %sum least-squares error

meanhappy = mean(data.happyscore);
re=sum((data.happyscore-meanhappy).^2);
happyr2 = 1-lse/re;