function [regressorScans,regressorNames] =  create_fMRIregressors_otherRPEactPass(tau, data, R, initialISIandSetup)
% an extract of fit_happy_model_otherRPEactPass(mtx,happyscore,constant)
% without fitting

if ~exist('tau','var'); error('Do use Tau'); end
plotFlag = 0;

%needs all blocks to be the same size
decayvec = 2*tau.^[0:size(data.certainmtx,2)-1]; decayvec = decayvec(:); % added *2 because we work here with 60 not 30 trials
cert = data.certainmtx; ev = data.evmtx; rpe = data.rpemtx; dec = decayvec;

theychose = data.theychosemtx; youchose = data.youchosemtx; cert = data.certainmtx; ev = data.evmtx; rpe = data.rpemtx; otherrpe = data.otherrpemtx; dec = decayvec;
otherrpeactive = youchose.*data.otherrpemtx;
otherrpepassive = theychose.*data.otherrpemtx;

regressorsTrials(:,1) = cert*dec;             regressorNames{1} = 'CR_withDecay';    % known from moment of decision
regressorsTrials(:,2) = ev*dec;               regressorNames{2} = 'EV_withDecay';     % known from moment of decision
regressorsTrials(:,3) = rpe*dec;              regressorNames{3} = 'selfRPE_withDecay';   % known from moment of outcome
regressorsTrials(:,4) = otherrpeactive*dec;   regressorNames{4} = 'otherRPEactive_withDecay'; % known from moment of outcome
regressorsTrials(:,5) = otherrpepassive*dec;  regressorNames{5} = 'otherRPEpassive_withDecay'; % known from moment of outcome
regressorsTrials = double(regressorsTrials);

if iscell(R); R = R{1}; end

resolution = 10;
TR = 2.5;
Ntr = 60; % N trials per run; that's the number a run is SUPPOSED to have

for s = 1:2 % separate by session
  NtrThisRun = Ntr;
  idx = [1:Ntr]+Ntr*(s-1);
  if idx(end) > length(R) % some trials missing
    idx = Ntr*(s-1)+1:length(R); % then take all remaining trials
    NtrThisRun = length(idx);
  end
  
  if ~isnan(initialISIandSetup(s))
    % get onsets
    startTime = R(idx(1)).optionsShown;
    lastOnset = R(idx(end)).optionsShown-startTime;
    fprintf('Last onset: %.0f ',lastOnset);
    onsetsOptions =   round(([R(idx).optionsShown] - startTime + initialISIandSetup(s))/TR)*resolution;
    %   onsetsChoice =    round(([R(idx).choiceShown] - startTime + initialISIandSetup(s))/TR)*resolution;
    outcomeTime =     round(([R(idx).outcomeShown] - startTime + initialISIandSetup(s))/TR)*resolution;
    RT =              round(nansum([[R(idx).acceptDecision]-[R(idx).optionsShown];[R(idx).RT]/1000]));
    decisionTime =    onsetsOptions + RT;
    
    % determine size of matrix to fill with regressorsTrials
    Nbins = round(lastOnset*1.05/TR*resolution);
    regressorsSubTR{s} = zeros(Nbins,5);
    
    % fill regressorsSubTR with regressorsTrials
    regressorsSubTR{s}(decisionTime(1:NtrThisRun),1:2) = regressorsTrials(idx,1:2);
    regressorsSubTR{s}(outcomeTime(1:NtrThisRun),3:5) = regressorsTrials(idx,3:5);
    
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
      keyboard
    end
  else
    regressorScans{s} = NaN;
  end
end
