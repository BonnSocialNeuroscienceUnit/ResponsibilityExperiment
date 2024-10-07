function [regressorScans,regressorNames] =  create_fMRIregressors(tau, data, R, initialISIandSetup)
% an extract of fit_happy_model_nogain_guiltenvy(mtx,happyscore,constant)
% without fitting

if ~exist('tau','var'); load('js_temp.mat'); end
if ~exist('tau','var'); initialISIandSetup = [3 3]; end
plotFlag = 0;

%needs all blocks to be the same size
decayvec = 2*tau.^[0:size(data.certainmtx,2)-1]; decayvec = decayvec(:); % added *2 because we work here with 60 not 30 trials
cert = data.certainmtx; ev = data.evmtx; rpe = data.rpemtx; dec = decayvec;
%ywtw = data.notguiltmtx;
ywtl = data.guiltmtx; %double(data.guiltmtx>0); %or leave as parametric - how much they got less than you
yltw = data.envymtx; %double(data.envymtx>0); %or leave as parametric - how much they got more than you
%yltl = data.notenvymtx;
regressorsTrials(:,1) = cert*dec;         regressorNames{1} = 'certain_rewards_withDecay';    % known from moment of decision
regressorsTrials(:,2) = ev*dec;           regressorNames{2} = 'expected_value_withDecay';     % known from moment of decision
regressorsTrials(:,3) = rpe*dec;          regressorNames{3} = 'prediction_error_withDecay';   % known from moment of outcome
regressorsTrials(:,4) = ywtl*dec;         regressorNames{4} = 'partner_got_x_less_withDecay'; % known from moment of outcome
regressorsTrials(:,5) = yltw*dec;         regressorNames{5} = 'partner_got_x_more_withDecay'; % known from moment of outcome
regressorsTrials = double(regressorsTrials);

if iscell(R); R = R{1}; end

resolution = 10;
TR = 2.5;
Ntr = 60; % N trials per run

for s = 1:2 % separate by session
  idx = [1:Ntr]+Ntr*(s-1); % idx of trials for this run
  
  if ~isnan(initialISIandSetup(s))
    % get onsets
    startTime = R(idx(1)).optionsShown;
    lastOnset = R(idx(end)).optionsShown-startTime;
    fprintf('Last onset in s: %.0f\n',lastOnset);
    onsetsOptions =   round(([R(idx).optionsShown] - startTime + initialISIandSetup(s))/TR)*resolution;
    %   onsetsChoice =    round(([R(idx).choiceShown] - startTime + initialISIandSetup(s))/TR)*resolution;
    outcomeTime =     round(([R(idx).outcomeShown] - startTime + initialISIandSetup(s))/TR)*resolution;
    RT =              round(nansum([[R(idx).acceptDecision]-[R(idx).optionsShown];[R(idx).RT]/1000]));
    decisionTime =    onsetsOptions + RT;
    
    % determine size of matrix to fill with regressorsTrials
    Nbins = round(lastOnset*1.05/TR*resolution);
    regressorsSubTR{s} = zeros(Nbins,5);
    
    % fill regressorsSubTR with regressorsTrials
    regressorsSubTR{s}(decisionTime(1:Ntr),1:2) = regressorsTrials(idx,1:2);
    regressorsSubTR{s}(outcomeTime(1:Ntr),3:5) = regressorsTrials(idx,3:5);
    
    % for each regressor: convolve with HRF
    for c = 1:size(regressorsSubTR{s},2)
      temp = conv(regressorsSubTR{s}(:,c),spm_hrf(TR/resolution)); % convolve at 10th TR resolution
      temp = temp(1:resolution:Nbins,:); % downsample to TR resolution
      regressorScans{s}(1:length(temp),c) = temp;
    end
    if plotFlag
      figure('position',[77 587 1113 702])
      subplot(3,1,1)
      plot(regressorsTrials(idx,:))
      xlabel('Trials')
      title(R(idx(1)).subjName,'interpreter','none')
      subplot(3,1,2)
      plot(regressorsSubTR{s})
      xlabel('Time in 10th of TRs')
      subplot(3,1,3)
      plot(regressorScans{s})
      xlabel('Time in TR')
      legend(regressorNames,'interpreter','none')
    end
  else
    regressorScans{s} = NaN;
  end
end
