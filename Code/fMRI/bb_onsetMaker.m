function [ons,dur,tr_names,values,missing_onsets] = onsetMaker(fname,startISI)
% Analysis script for Responsibility_fMRI
% rearranges the onsets of the important events of the fMRI experiment into
% a cell array according to their distinctive trial types
% Number of Trial types: 19

figFlag = 0;

load(fname);

% define basic variables to find trial types
cond = [results(:).condition]; % 1: social active; 2: social passive; 3: non-social
chRisky = [results(:).play];  chRisky(chRisky==0) = NaN;
chSafe = 1-[results(:).play];  chSafe(chSafe==0) = NaN;
CR = [results(:).safeOption] .* chSafe; % subject is certain to get this as the safe option was chosen
EV = mean(reshape([results(:).riskyOption],2,[])) .* chRisky; % average gain from risky options, which subject chose
PEself = EV-[results.amountThisTrialSubject]; % prediction error from real outcome of risky option
PEother = EV-[results.amountThisTrialPartner]; % prediction error from real outcome of risky option
outcomeDif = [results(:).amountThisTrialSubject]-[results(:).amountThisTrialPartner]; % the difference in the amount of money obtained by subject and partner in a given trial

% set reference time if none available
offset = results(1).optionsShown;
if ~exist('startISI','var')
  startISI = 0;
end

% all kinds of durations and onsets
RT =     [results(:).RT]/1000; RT(isnan(RT)) = 3; % response time, for opponent this was always 3s
onsOpt = [results(:).optionsShown] - offset + startISI; % time when options were shown
onsDec = [[results(:).optionsShown] + RT] - offset + startISI; % time when decision was taken or shown
onsOut = [results(:).outcomeShown] - offset + startISI; % time when the outcome was shown
onsHap = [results(:).ratingShown] - offset + startISI; % time when the happiness rating was given

% trial types: names, onsets, durations and parametric modulators (values)
% for options, onsets are the time at which the participant was presented the options
i = 1;  tr_names{i} = 'options_social';              idx = cond==1;     ons{i} = onsOpt(idx);   values{i} = CR(idx);
i = 2;  tr_names{i} = 'options_partner';             idx = cond==2;     ons{i} = onsOpt(idx);   values{i} = CR(idx);
i = 3;  tr_names{i} = 'options_solo';                idx = cond==3;     ons{i} = onsOpt(idx);   values{i} = CR(idx);
% for decisions, onsets are the time at which participant took the decision
i = 4;  tr_names{i} = 'decision_safe_social';        idx = ~isnan(chSafe) & cond==1;     ons{i} = onsDec(idx);   values{i} = CR(idx);
i = 5;  tr_names{i} = 'decision_safe_partner';       idx = ~isnan(chSafe) & cond==2;     ons{i} = onsDec(idx);   values{i} = CR(idx);
i = 6;  tr_names{i} = 'decision_safe_solo';          idx = ~isnan(chSafe) & cond==3;     ons{i} = onsDec(idx);   values{i} = CR(idx);
i = 7;  tr_names{i} = 'decision_risky_social';       idx = ~isnan(chRisky) & cond==1;    ons{i} = onsDec(idx);   values{i} = EV(idx);
i = 8;  tr_names{i} = 'decision_risky_partner';      idx = ~isnan(chRisky) & cond==2;    ons{i} = onsDec(idx);   values{i} = EV(idx);
i = 9;  tr_names{i} = 'decision_risky_solo';         idx = ~isnan(chRisky) & cond==3;    ons{i} = onsDec(idx);   values{i} = EV(idx);
% for outcomes, onsets are the time at which the participant was presented with the outcome
i = 10; tr_names{i} = 'outcome_safe_social';         idx = ~isnan(chSafe) & cond==1;     ons{i} = onsOut(idx);   values{i} = CR(idx);
i = 11; tr_names{i} = 'outcome_safe_partner';        idx = ~isnan(chSafe) & cond==2;     ons{i} = onsOut(idx);   values{i} = CR(idx);
i = 12; tr_names{i} = 'outcome_safe_solo';           idx = ~isnan(chSafe) & cond==3;     ons{i} = onsOut(idx);   values{i} = CR(idx);
i = 13; tr_names{i} = 'outcome_other_pos_social';    idx = ~isnan(chRisky) & cond==1 & PEother>0; ons{i} = onsOut(idx);   values{i} = zeros(1,length(ons{i}));
i = 14; tr_names{i} = 'outcome_other_pos_partner';   idx = ~isnan(chRisky) & cond==2 & PEother>0; ons{i} = onsOut(idx);   values{i} = zeros(1,length(ons{i}));
i = 15; tr_names{i} = 'outcome_self_pos_solo';       idx = ~isnan(chRisky) & cond==3 & PEself>0;  ons{i} = onsOut(idx);   values{i} = zeros(1,length(ons{i}));
i = 16; tr_names{i} = 'outcome_other_neg_social';    idx = ~isnan(chRisky) & cond==1 & PEother<0; ons{i} = onsOut(idx);   values{i} = zeros(1,length(ons{i}));
i = 17; tr_names{i} = 'outcome_other_neg_partner';   idx = ~isnan(chRisky) & cond==2 & PEother<0; ons{i} = onsOut(idx);   values{i} = zeros(1,length(ons{i}));
i = 18; tr_names{i} = 'outcome_self_neg_solo';       idx = ~isnan(chRisky) & cond==3 & PEself<0;  ons{i} = onsOut(idx);   values{i} = zeros(1,length(ons{i}));
% for happiness, onsets are the time at which the participant reported happiness - we don't have the onset of the happiness rating display
i = 19; tr_names{i} = 'happiness_rating';            idx = ~isnan([results(:).happiness]);   ons{i} = onsHap(idx);  values{i} = [results(idx).happiness];

% durations are all 0 for now.
for c = 1:length(tr_names)
  dur{c} = zeros(1,length(ons{c}));
end

% remove NaNs, round to 6 positions after comma
for t = 1:length(ons)
  ons{t} = ons{t}(~isnan(ons{t}));
  ons{t} = round(ons{t},6);
  dur{t} = dur{t}(~isnan(dur{t}));
end

% in case no onset available, replace with dummy very large onset
missing_onsets = zeros(1,length(ons));
for t = 1:length(ons)
  if isempty(ons{t})
    missing_onsets(t) = 1;
    ons{t} = max([ons{:}])-10*rand;
    dur{t} = 0;
    values{t} = 0;
  end
end
% disp(cellfun(@length,ons))

% to check: onsets should increase linearly, and as they are in seconds, if
% divided by 60 the max value should be around 20 as this was roughly the
% duration of the scan
if figFlag; figure; plot(sort([ons{:}]/60),'.'); ylabel('Time in min.'); xlabel('Onsets of all types');
  figure
  for c=1:length(ons)
    subplot(length(ons),1,c)
    stem(ons{c},ones(1,length(ons{c})))
  end
end
