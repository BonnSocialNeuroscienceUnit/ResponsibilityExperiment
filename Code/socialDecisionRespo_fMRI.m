function results = socialDecisionRespo_fMRI
% Code to run the fMRI experiment of the Responsibility Experiment. All
% instructions in German.
%
% J Schultz, M Gädeke and T Willems, University of Bonn 2017-2024
% Based on the PsychtoolBox 3's OldNewRecogExp.m (http://psychtoolbox.org)
%
% Builds on the experiment in Rutledge et al. https://doi.org/10.1038/ncomms11825.

fMRI = 1; % IsOSX = 1;

% Check for screen resulution, FHD doesn't work with 3T Scanner Monitor
if fMRI
    if get(0,'ScreenSize') ~= [1 1 1600 900]
        error('adjust screen resolution to 1600x900')
    end
end

% Get input from experimenters
prompt={'Please enter participant name:','Please enter partner name:'};
name='';
defaultanswer={'Respo-fMRI9999','Partner'};
answers=inputdlg(prompt,name,1,defaultanswer);
subjName=answers{1};
partnerName=answers{2};

% ---------- Preparation, buttons, stim size and timing -------------------
clc; rand('seed',sum(100*clock));

% Make sure keyboard mapping is the same on all supported operating systems
KbName('UnifyKeyNames');
escapeKeycode = KbName('ESCAPE');
spaceKeycode = KbName('SPACE');
triggerKeycode = KbName('s');
trigTimeManual = [];
% enterKeycode = KbName('RETURN'); enterKeycode = enterKeycode(1); % the first element is the RETURN key
if fMRI == 1
    rightRespKeycode = KbName('c'); leftRespKeycode = KbName('a');
    rightRespName = 'rechter Zeigefinger'; leftRespName = 'linker Zeigefinger';
    %     rightRespKeycode = KbName('8*'); leftRespKeycode = KbName('9(');
    %     rightRespName = 'Mittelfinger'; leftRespName = 'Zeigefinger';
    enterKeycode = KbName('d'); EnterKeyName = 'rechtem Daumen';
    buttonRefreshTimeB = 0.03;
else
    rightRespKeycode = KbName('RightArrow'); leftRespKeycode = KbName('LeftArrow');
    rightRespName = 'Rechter Pfeil'; leftRespName = 'Linker Pfeil';
    enterKeycode = KbName('SPACE'); enterKeycode = enterKeycode(1); % the first element is the RETURN key
    EnterKeyName = 'Leertaste';
    buttonRefreshTimeB = 0.05;
end

% ajust interval lenghts to fMRI/behavioural conditions
if fMRI == 1
    ISIs = 3:11;
else
    ISIs = 1:2;
end

scrSiz = get(0,'MonitorPositions');
fontSizeNormal = round(scrSiz(4)/25);    % adjust if screen size problem
fontSizeTitle = round(scrSiz(4)/40);
if ispc; lineHeight = 1.5; else; lineHeight = 1; end

results = struct; resEnd = [];

% ---------- Data files ---------------------------------------------------
% Define filenames of input files and result file:
if ~exist('data','dir')
    mkdir('data')
end

fs = filesep;

datafilename1 = sprintf('data%ssResponsibility_fMRI_%s_%s.mat',fs,subjName,datetime_js); % to save output

% ---------- Experiment setup ---------------------------------------------
% from Rutledge et al, PNAS 2014, SuppMat
NgainTrials = 20;
GainSafeGains = [20 30 40 50 60]/2; % divided to lower highest possible outcome
GainRiskyGainMultipliers = [1.68 1.82 2 2.22 2.48 2.8 3.16 3.6 4.2 5];
GainRiskyLosses = 0;

NmixedTrials = 20;
MixedSafeGains = 0;
MixedRiskyGains = [30 50 80 110 150]/2; % divided to lower highest possible outcome
MixedRiskyLossMultipliers = [0.2 0.34 0.5 0.64 0.77 0.89 1 1.1 1.35 2];

Nlosstrials = 20;
LossSafeGains = -GainSafeGains;
LossRiskyGains = 0;
LossRiskyLossMultipliers = GainRiskyGainMultipliers;

Ntrials = NgainTrials + NmixedTrials + NgainTrials;

startTime = [];

probabilites = round(gampdf(ISIs, 6,1)*100); %shape and scale paramaters from Rutledge et. al
ISI=[]; for i=1:length(ISIs); ISI=[ISI repmat(ISIs(i),1,probabilites(i))]; end

durSubjectDecides = 5;
durPartnerDecides = 3;
durShowChoiceMade = 2;
durPlayToOutcome = 1;
durOutcome = 2.5;

trialTypeGainMixedOrLoss = shuffle([ones(1,NgainTrials) 2*ones(1,NmixedTrials) 3*ones(1,Nlosstrials)]); % 1 is gain, 2 is mixed, 3 is loss
trialTypeCondition = shuffle([ones(1,length(trialTypeGainMixedOrLoss)/3) ...
    2*ones(1,length(trialTypeGainMixedOrLoss)/3) 3*ones(1,length(trialTypeGainMixedOrLoss)/3)]); % 1 = socially active trial, 2 = socially passive

% Instructions
message = sprintf([...
    'In diesem Experiment geht es um Entscheidungen, bei denen man Geld gewinnen oder verlieren kann.\n',...
    'In jeder Runde gibt es eine sichere und eine riskante Option.\n',...
    '\n',...
    'Sie und %s entscheiden abwechselnd fuer beide Spieler.\n',...
    '\n',...
    'Wenn Sie spielen, muessen Sie sich zwischen der sicheren und der riskanten Option entscheiden.\n',...
    'Waehlen Sie die sichere Option, bekommen Sie und ihr Partner den dazugehoerigen Betrag.\n',...
    'Waehlen Sie die riskante Option, bekommen beide mit gleicher Wahrscheinlichkeit einen der moeglichen Betraege.\n',...
    '\n',...
    'Wenn %s spielt, entscheidet er/sie fuer beide.\n',...
    'Falls die riskante Option gewaehlt wird, kann auch hier jeder separat den\n'...
    'besseren oder den schlechteren Betrag bekommen.\n',...
    '\n',...
    'Der Gewinn aller Spielrunden wird zusammengezaehlt und Ihnen beiden am Ende ausbezahlt.\n',...
    '\n',...
    'Ausserdem werden Sie ab und zu gefragt wie zufrieden Sie gerade sind.\n',...
    '\n',...
    'Bitte antworten Sie mit den Tasten "%s" und "%s".\n',...
    'Beginnen Sie das Experiment mit einem beliebigen Tastendruck.'],...
    partnerName,partnerName,leftRespName,rightRespName);
%    'Sie haben %d Sekunden um sich zu entscheiden.\n',...

% ---------- start PsychToolbox -----------------------------------------
try
    % prepare screen, hide cursor
    screens = Screen('Screens');
    screenNumber = max(screens);
    HideCursor;
    gray = GrayIndex(screenNumber);
    Screen('Preference', 'SkipSyncTests', 1)
    [w, wRect] = Screen('OpenWindow',screenNumber, gray);
    %   Screen('BlendFunction',w)
    Screen('TextSize', w, fontSizeTitle);
    offset = wRect(3)/5; % Determine left-right offset of stimuli based on screen size
    %   priorityLevel = MaxPriority(w);  Priority(priorityLevel);   % Set priority for script execution to realtime priority:
    % ListenChar(2) % Switch off character listening:
    
    messagePosTop = .3;        % the multiplication factor determining the position of the messages at the top of the screen
    messagePosBottom = 1.5;     % the multiplication factor determining the position of the messages at the bottom of the screen
    scrMidH = wRect(3)/2;
    scrMidV = wRect(4)/2;
    center = [scrMidH scrMidV];
    white = [255 255 255];
    
    % ----------- Setup variables slider ----------------------------------
    sliderFlankTextLeftB = 'gar nicht  |';
    sliderFlankTextRightB = '|  sehr';
    valenceMinMax = [1 51]; % so -1 and x2 -> ratings between 0 and 100
    sliderStepSizeB = 2;
    sliderIniB = repmat(' ',1,valenceMinMax(2));
    cursorB = 'O';
    message2{1} = 'Wie zufrieden sind Sie?';
    message2{2} = ['\nMarker bewegen mit Zeigefingern,\n bestaetigen mit ' EnterKeyName];
    sliderBoundsB = Screen('TextBounds',w,[sliderFlankTextLeftB sliderIniB sliderFlankTextRightB]);
    sliderBoundsB = sliderBoundsB(3:4);
    
    % ---------- All ready, show instructions -----------------------------
    % Show instructions to subject, centered and in white, wait for buttonpress to start:
    DrawFormattedText(w, message, 'center', 'center', white, [], [], [], lineHeight);  Screen('Flip', w);
    KeyIsDown = 0;  while ~KeyIsDown;    [KeyIsDown, startExpt] = KbCheck;  end
    Screen('Flip', w);
    Screen('TextSize', w, fontSizeNormal);
    
    % ---------- Scanner trigger --------------------------------------------
    KeyIsDown = 0; KeyCode = zeros(1,256);
    while ~KeyIsDown; [KeyIsDown, ~, KeyCode]=KbCheck;  end
    if ~fMRI  % Wait for button press to start:
        vbl = Screen('Flip', w);
    else % for fMRI, waits for button press to show Ready to scan screen
        vbl = Screen('Flip', w);
        % *************** FMRI TRIGGER WAITING HERE ***************
        % triggers and buttons
        
        startAtTriggerNo = 6;
        if ~IsOSX % if not a Mac then need serial port
            joker = '';
            sampleFreq = 120;
            baudRate = 115200;
            try; portSpec = FindSerialPort([], 1); catch; error('Please connect docking station'); end
            specialSettings = [];
            
            % Compute maximum input buffer size for 1 hour worth of triggers coming
            % in at a expected sampleFreq Hz with a size of at most 1 Bytes each:
            InputBufferSize = sampleFreq * 3600;
            
            % Assign an interbyte readtimeout which is either 15 seconds, or 10 times
            % the expected time between consecutive datapackets at the given sampleFreq
            % sampling frequency, whatever's higher. Could go higher or lower than
            % this, but this seems a reasonable starter: Will give code and devices
            % time to start streaming, but will prevent script from hanging longer than
            % 15 seconds if something goes wrong with the connection:
            readTimeout = max(10 * 1/sampleFreq, 15);
            
            % HACK: Restrict maximum timeout to 21 seconds. This is needed on Macintosh
            % computers, because at least OS/X 10.4.11 seems to have a bug which can
            % cause the driver to hang when trying to stop at the end of a session if
            % the timeout value is set higher than 21 seconds!
            readTimeout = min(readTimeout, 21);
            portSettings = sprintf('%s %s BaudRate=%i InputBufferSize=%i Terminator=0 ReceiveTimeout=%f ReceiveLatency=0.0001',...
                joker, specialSettings, baudRate, InputBufferSize, readTimeout);
            
            % Open port portSpec with portSettings, return handle:
            DrawFormattedText(w, 'Starting serial port link setup...\n', 'center', 'center', WhiteIndex(w), [], [], [], lineHeight); %center(1)*.75, center(2)*.3, 255);
            vbl = Screen('Flip', w); % now ready for trigger
            
            try myport = IOPort('OpenSerialPort', portSpec, portSettings); catch, myport = -1; end
            
            if myport >= 0 % then we have a serial port, BUT DOES NOT MEAN WE HAVE A CABLE CONNECTED TO IT!!
                %    fprintf('Link online: Hit a key on keyboard to start trigger recording, after that hit any key to finish trigger collection.\n');
                DrawFormattedText(w, 'Serial port link online, press space to start collecting triggers, escape to abort\n', 'center', 'center', WhiteIndex(w), [], [], [], lineHeight); %center(1)*.75, center(2)*.3, 255);
                vbl = Screen('Flip', w); % now ready for trigger
                KeyCode = zeros(1,256);  KeyIsDown = 0;
                while ~KeyCode(escapeKeycode) & ~KeyCode(spaceKeycode); [~, ~, KeyCode] = KbCheck;  end
                if KeyCode(escapeKeycode);  sca; ShowCursor; disp('aborted'); try IOPort('Close', myport); end; if debug; keyboard; else; return; end; end
                if KeyCode(spaceKeycode); end % now starts reading port
                IOPort('Purge', myport) % deletes triggers received since opening the port
                
                DrawFormattedText(w, 'Waiting for triggers; start fMRI sequence.\n', 'center', 'center', WhiteIndex(w), [], [], [], lineHeight); %center(1)*.75, center(2)*.3, 255);
                vbl = Screen('Flip', w); % now ready for trigger
                HideCursor
                
                % Start asynchronous background trigger collection and timestamping. Use
                % blocking mode for reading data -- easier on the system:
                asyncSetup = sprintf('%s BlockingBackgroundRead=1 StartBackgroundRead=1', joker);
                IOPort('ConfigureSerialPort', myport, asyncSetup);
                
                % Trigger reception started: From now on, the driver will read data from
                % the serial port, byte by byte. Whenever 1 Byte has been received from the
                % trigger mechanism, the driver assumes a trigger pulse is complete.
                % Reception of the first byte of a new packet is timestamped in GetSecs()
                % time, and the later IOPort('Read') calls will return those timestamps...
                
                nTriggersB4start = startAtTriggerNo; trigTime = [];
                trigStart = GetSecs;
                % tic
                trigStart = GetSecs;
                % tic
                while nTriggersB4start % there still are n triggers to wait for
                    % Wait blocking for a new data packet of 1 trigger byte from
                    % the serial port, then return the packet data as uint8's plus the
                    % GetSecs receive timestamp 'treceived' of the start of each packet:
                    [~, trigTime(end+1)] = IOPort('Read', myport, 1, 1);
                    trigTimeManual(end+1) = GetSecs;
                    nTriggersB4start = nTriggersB4start - 1/2;
                    DrawFormattedText(w, num2str(floor(nTriggersB4start)), 'center', 'center', WhiteIndex(w)); %center(1)*.75, center(2)*.3, 255);
                    Screen('Flip', w);
                end % no more triggers to wait for
                %    fprintf('%d triggers detected, starting experiment\n', count);
            else % no serial port
                disp('no serial port found! Try taking serial port cable out and putting it back in, restart Matlab')
                sca; ListenChar; keyboard
            end
        else % if it's a Mac
            DrawFormattedText(w, 'Waiting for scanner,\n press escape to abort, spacebar to start collecting triggers.\n', 'center', 'center', WhiteIndex(w), [], [], [], lineHeight); %center(1)*.75, center(2)*.3, 255);
            vbl = Screen('Flip', w); % now ready for trigger
            pause(.5)
            KeyCode = zeros(1,256); %KeyIsDown = 0;
            while ~KeyCode(escapeKeycode) && ~KeyCode(spaceKeycode); [~, ~, KeyCode] = KbCheck;  end
            if KeyCode(escapeKeycode);  sca; ShowCursor; disp('aborted'); if debug; sca; ListenChar; keyboard; else; return; end; end
            
            DrawFormattedText(w, 'Waiting for triggers; start fMRI sequence.\n', 'center', 'center', WhiteIndex(w), [], [], [], lineHeight); %center(1)*.75, center(2)*.3, 255);
            vbl = Screen('Flip', w); % now ready for trigger
            HideCursor
            pause(.5)
            
            nTriggersB4start = startAtTriggerNo; trigTime = []; trigTimeManual = [];  KeyCode = zeros(1,256);  KeyIsDown = 0;
            trigStart = GetSecs;
            % tic
            while nTriggersB4start % there still are n triggers to wait for
                [~, ~, KeyCode] = KbCheck;
                if KeyCode(escapeKeycode); sca; ShowCursor; disp('aborted'); if debug; sca; ListenChar; keyboard; else; return; end; end
                if KeyCode(triggerKeycode)
                    trigTimeManual(end+1) = GetSecs;
                    nTriggersB4start = nTriggersB4start - 1;
                    pause(.2)
                end % got input from scanner
                DrawFormattedText(w, num2str(floor(nTriggersB4start)), 'center', 'center', WhiteIndex(w));
                Screen('Flip', w);
                KeyCode = zeros(1,256);
            end % no more triggers to wait for
        end
    end % fMRI or not
    startTime = Screen('Flip', w);
    
    % ---------- initialize money at 0 at start -----------------------------
    moneyS = 0;
    moneyP = 0;
    
    % ---------- Loop through trials ----------------------------------------
    for trial = 1:length(trialTypeGainMixedOrLoss)
        
        % safe default value if partner isn't active
        acceptDecision = NaN;
        
        % wait a bit between trials
        WaitSecs(randsample(ISI,1)); % draws 1 ISI with set probabilities, mean ISI is 5.96 sec
        
        % ---- determine if social or not -------
        Condition = trialTypeCondition(trial); % social trial: 2 gambles; non-social trial: (only 1 gamble) now socially active trial!
        
        % ---- determine if gamble is on left or right ------------------------
        rightLeft = (round(rand)*2)-1;
        
        % ---- determine if gain, mixed or loss, compute safe and risky amounts -------
        clear safe risky
        switch trialTypeGainMixedOrLoss(trial)
            case 1 % gain trial
                safe = sampleFrom(1,GainSafeGains);
                risky = [safe*sampleFrom(1,GainRiskyGainMultipliers) sampleFrom(1,GainRiskyLosses)]; % gain and loss options
            case 2 % mixed trial
                safe = sampleFrom(1,MixedSafeGains);
                gambleGain = sampleFrom(1,MixedRiskyGains);
                risky = [gambleGain -gambleGain*sampleFrom(1,MixedRiskyLossMultipliers)]; % gain and loss options
            case 3 % loss trial
                safe = sampleFrom(1,LossSafeGains);
                risky = [sampleFrom(1,LossRiskyGains) safe*sampleFrom(1,LossRiskyLossMultipliers)]; % gain and loss options
        end
        % round the values to the nearest cent value
        safe = round(safe);
        risky = round(risky);
        
        % ------- show the options --------------------------------------------
        gambleStr = sprintf('Spiel: %d oder %d cent',risky(1),risky(2));
        safeStr = sprintf('Sicher: %d cent',safe);
        if rightLeft < 0 % then gamble on the left
            optionStr = [gambleStr '         ' safeStr];
        else % then gamble on the right
            optionStr = [safeStr '         ' gambleStr];
        end
        if Condition == 1 % socially active condition
            Screen('TextSize', w, fontSizeNormal*2)
            strName = ['Sie'];
            DrawFormattedText(w,strName, 'center', scrMidV/2, white, [], [], [], lineHeight*2); % show text outcome
            Screen('TextSize', w, fontSizeNormal);
            str = ['\nkoennen fuer sich und ' partnerName ' waehlen zwischen:\n' optionStr];
            if rightLeft < 0
                DrawFormattedText(w, sprintf('%s           %s',leftRespName,rightRespName), 'center', scrMidV*messagePosBottom, white);
                %                 DrawFormattedText(w, sprintf('Spielen (%s) oder nicht (%s)?',leftRespName,rightRespName), 'center', scrMidV*messagePosBottom, white);
            else
                DrawFormattedText(w, sprintf('%s           %s',leftRespName,rightRespName), 'center', scrMidV*messagePosBottom, white);
                %                 DrawFormattedText(w, sprintf('Nicht spielen (%s) oder spielen (%s)?',leftRespName,rightRespName), 'center', scrMidV*messagePosBottom, white);
            end
        elseif Condition == 2 % socially passive condition
            Screen('TextSize', w, fontSizeNormal*2)
            DrawFormattedText(w, partnerName, 'center', scrMidV/2, white, [], [], [], lineHeight*2); % show text outcome
            Screen('TextSize', w, fontSizeNormal);
            str = ['\n kann fuer beide waehlen zwischen:\n' optionStr];
        else
            Screen('TextSize', w, fontSizeNormal*2)
            strName = ['Sie'];
            DrawFormattedText(w,strName, 'center', scrMidV/2, white, [], [], [], lineHeight*2); % show text outcome
            Screen('TextSize', w, fontSizeNormal);
            str = ['\nkoennen NUR fuer sich waehlen zwischen:\n' optionStr];
            if rightLeft < 0
                %                 DrawFormattedText(w, sprintf('Spielen (%s) oder nicht (%s)?',leftRespName,rightRespName), 'center', scrMidV*messagePosBottom, white);
                DrawFormattedText(w, sprintf('%s           %s',leftRespName,rightRespName), 'center', scrMidV*messagePosBottom, white);
            else
                %                 DrawFormattedText(w, sprintf('Nicht spielen (%s) oder spielen (%s)?',leftRespName,rightRespName), 'center', scrMidV*messagePosBottom, white);
                DrawFormattedText(w, sprintf('%s           %s',leftRespName,rightRespName), 'center', scrMidV*messagePosBottom, white);
            end
        end
        DrawFormattedText(w,str, 'center', scrMidV/2, white, [], [], [], lineHeight*2); % show text outcome
        optionsShown = Screen('Flip', w);
        
        % ------- read or determine choice -----------------------------------------
        play = 0; % initialise to -1
        if Condition == 1 || Condition == 3 % non-social condition: subject decides
            KeyIsDown = 0;
            while ~KeyIsDown % no time limit for decision
                %            while (GetSecs - optionsShown) <= durSubjectDecides && ~KeyIsDown
                [KeyIsDown, resptime, keyCode] = KbCheck;
                if KeyIsDown
                    if keyCode(leftRespKeycode) && rightLeft == -1 % then gamble was on right, and subject chose certain amount
                        play = 1;
                    elseif keyCode(leftRespKeycode) && rightLeft == 1 % then gamble was on left, and subject chose the gamble
                        play = 0;
                    elseif keyCode(rightRespKeycode) && rightLeft == -1 % then gamble was on right, and subject chose the gamble
                        play = 0;
                    elseif keyCode(rightRespKeycode) && rightLeft == 1 % then gamble was on left, and subject chose certain amount
                        play = 1;
                    elseif keyCode(escapeKeycode)
                        sca; ShowCursor; disp('aborted');
                        if fMRI == 1
                            try IOPort('Close', myport);
                                disp('port closed')
                            catch, disp('no port found')
                            end;
                        end
                        return
                    else % do nothing
                    end % choices
                end % key is down
            end % wait for response
            Screen('Flip', w);
            
            % collect RT
            RT = round(1000*(resptime-optionsShown));
        elseif Condition == 2 % social condition: partner (algorithm) decides
            KeyIsDown = 0; % a decision has been made (  been 1 before implementing button press by subject)
            RT = NaN;
            WaitSecs(durPartnerDecides); % simulate the time the partner took to decide
            Screen('TextSize', w, fontSizeNormal*2)
            DrawFormattedText(w,partnerName, 'center', scrMidV/2, white, [], [], [], lineHeight*2); % show text outcome
            Screen('TextSize', w, fontSizeNormal);
            messagePartnerDecision = 'Ihr Partner hat gewaehlt. \n Druecken Sie eine beliebige Taste um fortzufahren.';
            str = ['\n kann fuer beide waehlen zwischen:\n' optionStr];
            DrawFormattedText(w,str, 'center', scrMidV/2, white, [], [], [], lineHeight*2);
            DrawFormattedText(w,messagePartnerDecision, 'center', scrMidV*messagePosBottom, white, [], [], [], lineHeight*2); % show text outcome
            Screen('Flip',w);
            while ~KeyIsDown;    [KeyIsDown, ~, KeyCode]=KbCheck; end
            acceptDecision = GetSecs;
            if safe(1) < mean(risky(1,:)); play = 1; % if absolute value of safe option less than half the mean of the risky ones, play.
            else; play = 0; % otherwise take safe option.
            end
        end
        
        % wait a bit before showing the choice made
        WaitSecs(randsample(ISI,1));
        
        % ------- show choice made --------------------------------------------
        if KeyIsDown == 1 % a decision was made
            if play % show the risky option
                chosenStr = gambleStr;
            else % show the safe option
                chosenStr = safeStr;
            end
            if Condition == 1 % socially active condition
                str = ['Sie haben fuer beide gewaehlt:\n' chosenStr];
            elseif Condition == 2 % socially passive condition
                str = [partnerName ' hat fuer beide gewaehlt:\n' chosenStr];
            else % non social condition
                str = ['Sie haben NUR fuer sich selbst gewaehlt:\n' chosenStr];
            end
        else
            str = 'Sie haben nicht gewaehlt. Beide bekommen die schlechteste Option.';
            RT = NaN;
        end
        DrawFormattedText(w,str, 'center', scrMidV, white, [], [], [], lineHeight*2); % show text outcome
        choiceShown = Screen('Flip', w);
        
        %WaitSecs(randsample(ISI,1));
        WaitSecs(durShowChoiceMade);
        
        % ---- determine outcome ----------------------------------------------
        winS = NaN; winP = NaN; amountThisTrialP = NaN; amountThisTrialS = NaN;
        if KeyIsDown == 1 % a decision was made
            if play % chose risky option
                if Condition == 2 % need a win and lose possibility for both
                    winP = round(rand); % does partner win?
                    winS = round(rand); % does subject win?
                    if winP; amountThisTrialP = max(risky); % higher value obtained
                    else; amountThisTrialP = min(risky); % lower value obtained
                    end
                    if winS; amountThisTrialS = max(risky); % higher value obtained
                    else; amountThisTrialS = min(risky); % lower value obtained
                    end
                elseif Condition == 1 % win/lose possiblity only for subject
                    winS = round(rand); % does subject win?
                    winP = round(rand);  %before: winP = NaN;
                    if winP; amountThisTrialP = max(risky); % higher value obtained -- before: nothing
                    else; amountThisTrialP = min(risky); % lower value obtained  -- before: nothing
                    end
                    if winS; amountThisTrialS = max(risky); % higher value obtained
                    else; amountThisTrialS = min(risky); % lower value obtained
                    end
                else
                    winS = round(rand); % does subject win?
                    if winS; amountThisTrialS = max(risky); % higher value obtained
                    else; amountThisTrialS = min(risky); % lower value obtained
                    end
                end
                
            else % safe amount given to all players
                amountThisTrialS = safe;
                if Condition < 3;
                    amountThisTrialP = safe;
                end
            end % play or not
        else % no decision made: worst outcome given.
            amountThisTrialS = min([safe risky]);
            amountThisTrialP = min([safe risky]);
        end
        
        % -------- show outcome -----------------------------------------------
        % wait a bit before outcome
        WaitSecs(durPlayToOutcome);
        
        if Condition == 1 % non-social condition: only 1 outcome
            str = sprintf('Sie bekommen %d cent, \n%s bekommt %d cent',amountThisTrialS, partnerName, amountThisTrialP);
        elseif Condition == 2 % social condition: 2 outcomes
            str = sprintf('%s bekommt %d cent, \nSie bekommen %d cent.',partnerName,amountThisTrialP,amountThisTrialS);
        else
            str = sprintf('Sie bekommen %d cent, \n%s ist davon nicht betroffen',amountThisTrialS, partnerName);
        end
        DrawFormattedText(w,str, 'center', scrMidV, white, [], [], [], lineHeight*2); % show text outcome
        outcomeShown = Screen('Flip', w);
        
        WaitSecs(durOutcome);
        Screen('Flip', w);
        
        % -------- increment subject gains with amount ------------------------
        moneyS = moneyS + amountThisTrialS;
        moneyP = nansum([moneyP,amountThisTrialP]);
        
        % -------- get happiness rating ---------------------------------------
        getHappiness = 1-rem(trial,2);
        happiness = NaN;
        if getHappiness
            % prepare slider variables for response
            slider = sliderIniB; slider(round(length(slider)/2)) = cursorB;
            resp = round(length(slider)/2);
            str = [message2{1} '\n' sliderFlankTextLeftB slider sliderFlankTextRightB '\n' message2{2}];
            DrawFormattedText(w,str, 'center', scrMidV-scrMidV/2, white, [], [], [], lineHeight*2); % show text outcome
            ratingShown = Screen('flip', w);
            ratingGiven = 0; KeyIsDown = 0;
            while ~ratingGiven
                [KeyIsDown, resptime, keyCode] = KbCheck;
                if KeyIsDown
                    if keyCode(leftRespKeycode) & resp > 1
                        resp = resp - sliderStepSizeB;
                    elseif keyCode(rightRespKeycode) & resp < length(slider)
                        resp = resp + sliderStepSizeB;
                    elseif any(keyCode(enterKeycode))
                        ratingGiven = 1;
                    elseif keyCode(escapeKeycode)
                        sca; ShowCursor; disp('aborted');
                        try IOPort('Close', myport);
                            disp('port closed')
                        catch, disp('no port found')
                        end;
                        return
                    else % do nothing
                    end
                    if resp<1; resp=1; end
                    if resp>valenceMinMax(2); resp=valenceMinMax(2); end
                    slider = sliderIniB;
                    slider(resp) = cursorB;
                    str = [message2{1} '\n' sliderFlankTextLeftB slider sliderFlankTextRightB '\n' message2{2}];
                    DrawFormattedText(w,str, 'center', scrMidV-scrMidV/2, white, [], [], [], lineHeight*2); % show text outcome
                    Screen('flip', w);
                    WaitSecs(buttonRefreshTimeB)
                end % collect button press entry
            end % got a rating
            happiness = resp;
            ratingEnd = GetSecs;
        else
            ratingShown = NaN;
            ratingEnd = NaN;
        end % end happiness rating
        Screen('Flip', w);
        
        % ---------- Finish task, collect answers and variables -------------
        % put trial information into results structure:
        results(trial).safeOption = safe; % how much was the safe amount?
        results(trial).optionsShown = optionsShown;
        results(trial).acceptDecision = acceptDecision;
        results(trial).choiceShown = choiceShown;
        results(trial).outcomeShown = outcomeShown;
        results(trial).ratingShown = ratingShown;
        results(trial).riskyOption = risky; % how much was the safe amount?
        results(trial).amountThisTrialPartner = amountThisTrialP; % how much did the partner get?
        results(trial).amountThisTrialSubject = amountThisTrialS; % how much did the subject get?
        results(trial).moneyP = moneyP;% cumulative earnings of partner
        results(trial).moneyS = moneyS;% cumulative earnings of subject
        results(trial).play = play;         % did subject or partner want to play or not?
        results(trial).RT = RT;         % what was the response time?
        results(trial).subjectWon = winS;    % did subject win?
        results(trial).partnerWon = winP;    % did partner win?
        results(trial).happiness = happiness;    % happiness rating
        results(trial).RThappiness = (ratingEnd-ratingShown)*1000; % duration of happiness rating
        results(trial).condition = Condition; % 1 = socially active; 2 = socially passive; 3 = non-social
        
        % ---------- Save data ------------------------------------------------
        if ~strcmpi(subjName,'demo') % only if not demo
            save(datafilename1,'results','startTime')      % Write results structure to file:
        end
        
        Screen('Flip', w);
        
    end % for trial loop
    endExpt = GetSecs;
    
    %Report earnings to subject
    str = sprintf('Das Ergebnis ihres Spiels:\n Sie erhalten %d cent',moneyS);
    DrawFormattedText(w,str, 'center', scrMidV, white, [], [], [], lineHeight*2); % show text outcome
    Screen('flip', w);
    WaitSecs(3)
    while ~KeyIsDown;    [KeyIsDown, ~, KeyCode]=KbCheck; end
    
    str = sprintf('Das Ergebnis ihres Partners:\n %s erhaelt %d cent',partnerName, moneyP);
    DrawFormattedText(w,str, 'center', scrMidV, white, [], [], [], lineHeight*2); % show text outcome
    Screen('flip', w);
    WaitSecs(3)
    while ~KeyIsDown;    [KeyIsDown, ~, KeyCode]=KbCheck; end
    % ----------- Cleanup at end of experiment ------------------------------
    save(datafilename1,'results','startTime')      % Write results structure to file:
    sca;
    ShowCursor;
    %    ListenChar;
    try IOPort('Close', myport);
        disp('port closed')
    catch disp('no port found')
    end;
    
    fprintf('Partner erhaelt %d Cent', results(end).moneyP);
    fprintf('Proband erhaelt %d Cent',results(end).moneyS);
    fprintf('Duration: %.0f seconds\n',endExpt-startExpt)
catch errMsg
    % Do same cleanup as at the end of a regular session...
    save(datafilename1,'results','startTime')      % Write results structure to file:
    sca;
    ShowCursor;
    % ListenChar;
    if fMRI == 1
      try IOPort('Close', myport);
        disp('port closed')
      catch disp('no port found')
      end;
    end
    keyboard
    % Output the error message that describes the error:
    psychrethrow(psychlasterror);
end % try ... catch %
