function [result] = fit_happy_model_Omar_Social(mtx,happyscore,constant)

%assumes that each row of the cell array is a different subject

rptfit = 100; %repeat N times jittering the betas
tempf=fieldnames(mtx);
for n=1:length(tempf), eval(sprintf('mtx.%s=double(mtx.%s);',tempf{n},tempf{n})); end;
data = mtx;
data.happyscore = double(happyscore);

%result = data;
if isfield(mtx,'inx')
    inx = mtx.inx;
else
    inx = [0.5 0.5 0.5 0.005 10 0.5 0.5]; %dummies 
end
options = optimset('Display','off','MaxIter',100000,'TolFun',1e-10,'TolX',1e-10,...
    'DiffMaxChange',1e-2,'DiffMinChange',1e-4,'MaxFunEvals',10000,'LargeScale','off');
warning off; %display,iter to see outputs
lb = [-100 -100 -100 0.005 0.1  -100 -100];
ub = [100   100  100 0.99  100    100 100];

if exist('constant','var') && ~isnan(constant)
    inx = [inx mean(happyscore)];
    lb = [lb 0];
    ub = [ub 1];
    if isfield(mtx,'const') %use previous fit constant
        inx(end) = mtx.const;
    end;
end;

result.modelName = 'minEV+PE';
result.inx = double(inx);
dof = length(inx);
%result.options = options;
result.lb  = lb;
result.ub  = ub;
b=inx;
%b = fmincon(@model_param, inx, [],[],[],[],lb,ub,[], options, data);
result.b = b;
[lse, happypred, r2] = model_param(b, data);

for n=1:rptfit %jitter the betas
    b = fmincon(@model_param, double(b+0.1*randn(size(b))), [],[],[],[],lb,ub,[], options, data);
    [lse2, happypred, r2] = model_param(b, data);
    if lse2<lse
        result.b = b;
       [lse, happypred, r2] = model_param(b, data);
    end;
end;

result.happypred = happypred;
result.r2 = r2;
result.bic = length(happypred)*log(lse/length(happypred)) + length(b)*log(length(happypred));
result.aic = length(happypred)*log(lse/length(happypred)) + 2*length(b);
result.paramNames = {'Ch_diff','MinEV','PE','tau','Utl','Guilt','Envy'};
if exist('constant','var') && ~isnan(constant)
  result.paramNames = [result.paramNames, {'constant'}];
end
end

function [lse, happypred, happyr2] =  model_param(x, data)

a = x(1);   %accumulated gain or loss
b = x(2);   %guranteed win paprameter
c = x(3);   %expected win parameter
tau = x(4); %decay constant
utl= x(5);    %Prediction error parameter
try s1 = x(6);  %guilt
  s2 = x(7);  %envy
end
% % % d= x(6);
% % % e= x(7);
% % % f= x(8);
% % % const= x(7);
% if length(x)==5
%     const = x(5);
% else
%     const = 0; %if z-scored
% end;

%needs all blocks to be the same size
decayvec = tau.^[0:size(data.Agl_Abs,2)-1]; 
% % % decayvec = [decayvec(1)./2 decayvec(1)./2 decayvec(2:end-1)]; 
decayvec=decayvec(:);
% figure;plot(decayvec)
happypred = data.happyscore;
theychose = data.theychosemtx; youchose = data.youchosemtx; cert = data.certainmtx; ev = data.evmtx; rpe = data.rpemtx; dec = decayvec;
ywtl = data.guiltmtx; %double(data.guiltmtx>0); %or leave as parametric - how much they got less than you
yltw = data.envymtx; %double(data.envymtx>0); %or leave as parametric - how much they got more than you
% theychose = data.theychosemtx; youchose = data.youchosemtx; 
agl = data.Agl_Abs; minEV = data.MinEV; ch_diff = data.Ch_diff; 
dec = decayvec; 
PE = data.PE;
% % % CR= data.certainmtx;
% % % GBEV= data.rewardmtx;
% M1=(( minEV + agl) ./ agl);
% M1n=find(isnan(M1));
% M1Inf=find(M1==Inf);
% M1Inf2= find(M1==-Inf);
% invalid_M1={M1n,M1Inf,M1Inf2};
% for invalid_type = 1:size(invalid_M1,2)
% if ~isempty(invalid_M1{invalid_type})
% for fill_M1_0= 1:size(invalid_M1{invalid_type},1)
% M1(invalid_M1{1,invalid_type}(fill_M1_0,1))=0;
% end
% end
% end
% M2=((ch_diff)./(minEV)).*100;
% M2n=find(isnan(M2));
% M2Inf=find(M2==Inf);
% M2Inf2= find(M2==-Inf);
% invalid_M2={M2n,M2Inf,M2Inf2};
% for invalid_type = 1:size(invalid_M2,2)
% if ~isempty(invalid_M2{invalid_type})
% for fill_M2_0= 1:size(invalid_M2{invalid_type},1)
% M2(invalid_M2{1,invalid_type}(fill_M2_0,1))=0;
% end
% end
% end
% M3=((PE)./(minEV+ch_diff)).*100;
% M3_5=find(isnan(M3));
% M3Inf=find(M3==Inf);
% M3Inf2= find(M3==-Inf);
% invalid_M3={M3_5,M3Inf,M3Inf2};
% for invalid_type = 1:size(invalid_M3,2)
% if ~isempty(invalid_M3{invalid_type})
% for fill_M3_0= 1:size(invalid_M3{invalid_type},1)
% M3(invalid_M3{1,invalid_type}(fill_M3_0,1))=0;
% end
% end
% end
% %  happypred = a*agl*dec + b*M1*dec + c*M2*dec + d*M3*dec;
% dec2=zeros(size(minEV,2),1);
% dec2(1,1)=1;
% Utl_minEV=1./(1+exp(-2.3.*(minEV-0)));
% Utl_minEV_5=find(Utl_minEV==0.5);
% invalid_M3={Utl_minEV_5,M3Inf,M3Inf2};
% for Utl_minEV_i=1:size(Utl_minEV_5,1)
% Utl_minEV(Utl_minEV_5(Utl_minEV_i))=0;
% end

range=14;
slope=0;

%start Making your model terms
% max_ch=max(ch_diff,1);
utl2=(utl.^2);
term1x=(ch_diff.^2)*dec;
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
% min_minEV=min(minEV,1);
term2x=minEV*dec;
term2x2=zeros(size(term2x,1),1);
for res_loi=1:size(term2x2,1)
if term2x(res_loi,1)<0 
% term2x2(res_loi,1)=-log(-term2x(res_loi,1))- 2.*utl.*((range/2)+(exp(term2x(res_loi,1))/2));% Minus Utility function
% term2x2(res_loi,1)=real(term2x2(res_loi,1));
term2x2(res_loi,1)=b.*((range./(1+exp(-(3/4).*(utl).*(term2x(res_loi,1)))))-(range/2));
if term2x2(res_loi,1)< -(range/2)
term2x2(res_loi,1)= -(range/2);
elseif term2x2(res_loi,1)==-Inf
  term2x2(res_loi,1)= -(range/2);
end
elseif term2x(res_loi,1)>0
% term2x2(res_loi,1)=log(term2x(res_loi,1))+ utl.*((range/2)+(exp(-term2x(res_loi,1))/2)); % Positive Utility function
term2x2(res_loi,1)=b.*((range./(1+exp(-(3/12).*utl.*(term2x(res_loi,1)))))-(range/2));
if term2x2(res_loi,1)> (range/2)
    term2x2(res_loi,1)= (range/2);
end
end
end
% term3x=0;
% min_PE=min(PE,1);
% max_PE=max(PE,1);
% if -(min_PE)>max_PE
% tail_PE=min_PE;
% term3x=(3.*(1-(PE./tail_PE)))*dec;
% else
%   tail_PE=max_PE;  
% term3x=(3.*PE./tail_PE)*dec;
% end
term3x=PE*dec;
term3x2=zeros(size(term3x,1),1);
for res_loi=1:size(term3x2,1)
if term3x(res_loi,1)<0
% term3x2(res_loi,1)=-log(-term3x(res_loi,1))- 2.*utl.*((range/2)+(exp(term3x(res_loi,1))/2)); %Minus Utility function
% term3x2(res_loi,1)=real(term3x2(res_loi,1));
term3x2(res_loi,1)=c.*((range./(1+exp(-(3/4).*(utl2).*(term3x(res_loi,1)))))-(range/2));
if term3x2(res_loi,1)<= -(range/2)
term3x2(res_loi,1)= -(range/2);
elseif term3x2(res_loi,1)==-Inf
  term3x2(res_loi,1)= -(range/2);
end
elseif term3x(res_loi,1)>0
% term3x2(res_loi,1)=log(term3x(res_loi,1))+ utl.*((range/2)+(exp(-term3x(res_loi,1))/2)); %Utility function
term3x2(res_loi,1)=c.*((range./(1+exp(-(3/12).*utl2.*(term3x(res_loi,1)))))-(range/2));
if term3x2(res_loi,1)>= (range/2)
    term3x2(res_loi,1)= (range/2);
end
end
end
term4x=ywtl*dec;
term4x2=zeros(size(term4x,1),1);
for res_loi=1:size(term4x2,1)
if term4x(res_loi,1)<0
% term3x2(res_loi,1)=-log(-term3x(res_loi,1))- 2.*utl.*((range/2)+(exp(term3x(res_loi,1))/2)); %Minus Utility function
% term3x2(res_loi,1)=real(term3x2(res_loi,1));
term4x2(res_loi,1)=s1.*((range./(1+exp(-(3/4).*(utl2).*(term4x(res_loi,1)))))-(range/2));
if term4x2(res_loi,1)<= -(range/2)
term4x2(res_loi,1)= -(range/2);
elseif term4x2(res_loi,1)==-Inf
  term4x2(res_loi,1)= -(range/2);
end
elseif term4x(res_loi,1)>0
% term3x2(res_loi,1)=log(term3x(res_loi,1))+ utl.*((range/2)+(exp(-term3x(res_loi,1))/2)); %Utility function
term4x2(res_loi,1)=s1.*((range./(1+exp(-(3/12).*utl2.*(term4x(res_loi,1)))))-(range/2));
if term4x2(res_loi,1)>= (range/2)
    term4x2(res_loi,1)= (range/2);
end
end
end
term5x=yltw*dec;
term5x2=zeros(size(term5x,1),1);
for res_loi=1:size(term5x2,1)
if term5x(res_loi,1)<0
% term3x2(res_loi,1)=-log(-term3x(res_loi,1))- 2.*utl.*((range/2)+(exp(term3x(res_loi,1))/2)); %Minus Utility function
% term3x2(res_loi,1)=real(term3x2(res_loi,1));
term5x2(res_loi,1)=s2.*((range./(1+exp(-(3/4).*(utl2).*(term5x(res_loi,1)))))-(range/2));
if term5x2(res_loi,1)<= -(range/2)
term5x2(res_loi,1)= -(range/2);
elseif term5x2(res_loi,1)==-Inf
  term5x2(res_loi,1)= -(range/2);
end
elseif term5x(res_loi,1)>0
% term3x2(res_loi,1)=log(term3x(res_loi,1))+ utl.*((range/2)+(exp(-term3x(res_loi,1))/2)); %Utility function
term5x2(res_loi,1)=s2.*((range./(1+exp(-(3/12).*utl2.*(term5x(res_loi,1)))))-(range/2));
if term5x2(res_loi,1)>= (range/2)
    term5x2(res_loi,1)= (range/2);
end
end
end
% % % max_agl=max(agl,1);
% % % term4x=agl*dec;
% % % term4x2=zeros(size(term4x,1),1);
% % % for res_loi=1:size(term4x2,1)
% % % if term4x(res_loi,1)<0
% % % % term1x2(res_loi,1)=-log(-term1x(res_loi,1))- 2.*utl.*((range/2)+(exp(term1x(res_loi,1))/2));%Minus Utility function
% % % % term1x2(res_loi,1)=real(term1x2(res_loi,1));
% % %  term4x2(res_loi,1)=d.*((range./(1+exp(-(3/4).*(utl2).*(term4x(res_loi,1)))))-(range/2));
% % % if term4x2(res_loi,1)< -(range/2)
% % % term4x2(res_loi,1)= -(range/2);
% % % elseif term4x2(res_loi,1)==-Inf
% % %   term4x2(res_loi,1)= -(range/2);
% % % end
% % % elseif term4x(res_loi,1)>0
% % % % term1x2(res_loi,1)=log(term1x(res_loi,1))+ utl.*((range/2)+(exp(-term1x(res_loi,1))/2)); %Positive Utility function
% % %  term4x2(res_loi,1)=d.*((range./(1+exp(-(3/12).*(utl2).*(term4x(res_loi,1)))))-(range/2));
% % % if term4x2(res_loi,1)> (range/2)
% % %     term4x2(res_loi,1)= (range/2);
% % % end
% % % end
% % % end
% % % max_CR=max(CR,1);
% % % term5x=CR*dec;
% % % term5x2=zeros(size(term5x,1),1);
% % % for res_loi=1:size(term5x2,1)
% % % if term5x(res_loi,1)<0
% % % % term1x2(res_loi,1)=-log(-term1x(res_loi,1))- 2.*utl.*((range/2)+(exp(term1x(res_loi,1))/2));%Minus Utility function
% % % % term1x2(res_loi,1)=real(term1x2(res_loi,1));
% % %  term5x2(res_loi,1)=e.*((range./(1+exp(-(3/4).*(utl2).*(term5x(res_loi,1)))))-(range/2));
% % % if term5x2(res_loi,1)< -(range/2)
% % % term5x2(res_loi,1)= -(range/2);
% % % elseif term5x2(res_loi,1)==-Inf
% % %   term5x2(res_loi,1)= -(range/2);
% % % end
% % % elseif term5x(res_loi,1)>0
% % % % term1x2(res_loi,1)=log(term1x(res_loi,1))+ utl.*((range/2)+(exp(-term1x(res_loi,1))/2)); %Positive Utility function
% % %  term5x2(res_loi,1)=e.*((range./(1+exp(-(3/12).*(utl2).*(term5x(res_loi,1)))))-(range/2));
% % % if term5x2(res_loi,1)> (range/2)
% % %     term5x2(res_loi,1)= (range/2);
% % % end
% % % end
% % % end
% % % gbev_max=max(GBEV,1);
% % % term6x=GBEV*dec;
% % % term6x2=zeros(size(term6x,1),1);
% % % for res_loi=1:size(term6x2,1)
% % % if term6x(res_loi,1)<0
% % % % term1x2(res_loi,1)=-log(-term1x(res_loi,1))- 2.*utl.*((range/2)+(exp(term1x(res_loi,1))/2));%Minus Utility function
% % % % term1x2(res_loi,1)=real(term1x2(res_loi,1));
% % %  term6x2(res_loi,1)=f.*((range./(1+exp(-(3/4).*(utl2).*(term6x(res_loi,1)))))-(range/2));
% % % if term6x2(res_loi,1)< -(range/2)
% % % term6x2(res_loi,1)= -(range/2);
% % % elseif term6x2(res_loi,1)==-Inf
% % %   term6x2(res_loi,1)= -(range/2);
% % % end
% % % elseif term6x(res_loi,1)>0
% % % % term1x2(res_loi,1)=log(term1x(res_loi,1))+ utl.*((range/2)+(exp(-term1x(res_loi,1))/2)); %Positive Utility function
% % %  term6x2(res_loi,1)=f.*((range./(1+exp(-(3/12).*(utl2).*(term6x(res_loi,1)))))-(range/2));
% % % if term6x2(res_loi,1)> (range/2)
% % %     term6x2(res_loi,1)= (range/2);
% % % end
% % % end
% % % end
happypred=(term1x2+term2x2+term3x2+term4x2+term5x2)./5;
%happypred=happypredx*decx;
%b.*((-(ch_diff.^2)./1)+(range./(2+2.*(exp(-slope2.*ch_diff))))-(range/4))*dec
meanhappy = mean(data.happyscore);
lse = sum((data.happyscore-happypred).^2); %sum least-squares error (SSR)
% SSE= sum((happypred-meanhappy).^2); %SSE
% SST= sum((data.happyscore-meanhappy).^2); %SST
re=sum((data.happyscore-meanhappy).^2);
happyr2 = 1-(lse/re);
% lse = sum((happypred-mean(data.happyscore)).^2);  % Sum of squared error SSR
% TSS = sum((data.happyscore-mean(data.happyscore)).^2);     % Total sum of squares.
% happyr2 = SSE/SST;
end