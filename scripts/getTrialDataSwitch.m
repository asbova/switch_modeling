%% Code built to pull switch data from MedPC files
%% Austin Bruce, Tomas Lence, and Matt Weber

%% Changes to make - notes go here

%% Not sure what this is for - remnant of Ben's 
%change choice configuration list to 12 and 21 for left solid / right blink and left blink / right solid
%Add G as a variable and SET G(C(^trialCounter)) = 1 for left choices and SET G(C(^trialCounter)) = 2; for right pokes make sure everything is sealed
% S8, \ correction trial
%     0.01": 
%         IF (C(^goodTrial) >= 1)  [@LOGDATA;@IGNOREDATA] \ logging correct/incorrect/correction like this will make it easier to integrate with matlab (i.e., goodtrialsindex=L>0)
% 	   @LOGDATA:    SET L(C(^trialCounter)) = 2; ---> SX; \ log that this non-correction trial was correct
% 	   @IGNOREDATA: SET L(C(^trialCounter)) = 0; ---> SX; !!!\ log that this was a correction trial. 
% 
%         ADD C(^trialCounter); 
%         SET C(^CorrectionTrial) = 1 ---> S9 

%% following are data format from MEDPC output file
%{
From Bisection task

\   A       Probe trial selection array
\   B       Reward counter
\   D
\   E       Interval Record Array
\   F
\   G       center nosepoke in array
\   H       Trial Start time record array
\   I       Variable Interval arrary
\   K       List of random ITI range
\   L       Log good trial accuracy data
\   M       Correction Trial tracker

\   N	    RIGHT NP response time record array
\   O       RIGHT NP release time record array

\   P       LEFT NP Response Time (Beam Break) record array
\   Q       LEFT NP Response Time (Beam Release) record array
\   R	    Opto laser on (1) or off (2)
\   S       Trial End time record array
\   T       Probe trial record array
\   U       ITI record array
\   V	    Probe length
\   W       Pellet dispense time record array
\   X       milli second Timer for recoring events
\   Y       Reward zone out record array
\   Z       Reward zone in record array
\-----------------------------------------
\----------------------------------------------------------
%}

%% Starts to pull values from MedPC files
function TrialAnSt = getTrialDataSwitch(mpcParsed)
uniqueSbj = unique({mpcParsed.Subject});
ompcType = cellfun(@length, {mpcParsed.Q}) == 1 & cellfun(@length, {mpcParsed.Y}) == 0 & cellfun(@length, {mpcParsed.C}) == 0 & cellfun(@length, {mpcParsed.S}) > 0; % old MPC with trials
nmpcType = ~ompcType & min(cellfun(@length, {mpcParsed.H}), cellfun(@length, {mpcParsed.S}))> 0; % check if any trials
for i_animal = 1:size(uniqueSbj,2)
    sbjIdx = strcmp(uniqueSbj(i_animal),{mpcParsed.Subject}) & (ompcType | nmpcType);
    sbjLineIdx = find(sbjIdx);
    num_days = sum(sbjIdx);
    animals = char(uniqueSbj(i_animal));  
    if ~isempty(animals)
        for i_day = 1:num_days
            lineIdx =sbjLineIdx(i_day);
            
            if contains(mpcParsed(lineIdx).MSN,'18L6R')
                type = 1;
            elseif contains(mpcParsed(lineIdx).MSN,'6L18R')
                type = 0;
            else
                type = NaN;
            end

            timeITI = mpcParsed(lineIdx).U';
            trialStart = mpcParsed(lineIdx).H'; % 1
            trialEnd = mpcParsed(lineIdx).S'; % 2
            trialStart=trialStart(1:length(trialEnd));
            trialType = mpcParsed(lineIdx).T'; % the programmed duration for the trial
            trialType=trialType(1:length(trialEnd));
            reward = mpcParsed(lineIdx).W';
            opto = mpcParsed(lineIdx).R'; 
            opto = opto==1; % 1 means laser is on. 2 is off

            rpInLeft = mpcParsed(lineIdx).P'; % responding on left NP
            rpOutLeft = mpcParsed(lineIdx).Q';
            rpInRight = mpcParsed(lineIdx).N'; % responding on right NP
            rpOutRight = mpcParsed(lineIdx).O';
            rpInBack = mpcParsed(lineIdx).F'; % responding on right NP
            rpOutBack = mpcParsed(lineIdx).G'; % correct: being used for NP out and side choice...
            fzIn = mpcParsed(lineIdx).Z'; % enter food hopper
            fzOut = mpcParsed(lineIdx).Y'; % exit food hopper

            trialOutcomeInfo = mpcParsed(lineIdx).L'; % 4 = correct trial, 3 = incorrect trial, 2 = correct correction trial, 1 = incorrect correction trial

            trial = struct;
            trialNum = min(length(trialStart), length(trialEnd));
            for i_trial = 1:trialNum
                %% trial start/end time, duration, and ITI
                curTS = trialStart(i_trial); % H 
                curTE = trialEnd(i_trial);
                trial(i_trial).trialStart = curTS;
                trial(i_trial).trialEnd = curTE;
                trial(i_trial).trialDuration = curTE - curTS;
                trial(i_trial).ITI = timeITI(i_trial);

                %% real trial start time (when mouse iniated)
                x = rpInBack-curTS;
                trial(i_trial).initiationRT = min(x(x>=0)); % timestamp of first response after the trial was initiated
                curTS =  trial(i_trial).initiationRT +  curTS; % add the reaction time to trial start (when back light turned on), get the real trial start. 
                trial(i_trial).realTrialStart = curTS;

                %% programmed duration for the trial
                trial(i_trial).programmedDuration = trialType(i_trial);

                %% opto
                 if numel(opto)~=0; trial(i_trial).opto = opto(i_trial); end

                %% outcome info
                trial(i_trial).outcome = trialOutcomeInfo(i_trial);

                %% response time within trial
                trial(i_trial).leftRspTimeTrial = rpInLeft(rpInLeft>curTS & rpInLeft<=curTE) - curTS;
                trial(i_trial).leftRelTimeTrial = rpOutLeft(rpOutLeft>curTS & rpOutLeft<=curTE) - curTS;
                trial(i_trial).rightRspTimeTrial = rpInRight(rpInRight>curTS & rpInRight<=curTE) - curTS;
                trial(i_trial).rightRelTimeTrial = rpOutRight(rpOutRight>curTS & rpOutRight<=curTE) - curTS;

                %% Builds the Switch Departure and Switch Arrival Times based on behavioral training paradigm
                % 1. Start at short and then switch / 2. Start at short and then switch, back to short / 3. Multiple switches
                % 4. Start at long, move to short and then switch / 5. Start and stay at long / 6. Start and stay at short  
                if type == 0    % 0 indicates that a short latency trial is rewarded at the left nose poke 
                    if ~isempty(trial(i_trial).leftRelTimeTrial) && ~isempty(trial(i_trial).rightRspTimeTrial) 
                        trial(i_trial).SwitchArrival = min(trial(i_trial).rightRspTimeTrial(trial(i_trial).rightRspTimeTrial > min(trial(i_trial).leftRelTimeTrial)));
                            if ~isempty(trial(i_trial).SwitchArrival) 
                                trial(i_trial).SwitchDepart = max(trial(i_trial).leftRelTimeTrial(trial(i_trial).leftRelTimeTrial < trial(i_trial).SwitchArrival));
                            else trial(i_trial).SwitchArrival=[]; 
                            end
                    elseif isempty(trial(i_trial).leftRelTimeTrial) && isempty(trial(i_trial).rightRspTimeTrial)
                        trial(i_trial).SwitchArrival = [];
                        trial(i_trial).SwitchDepart = [];
                    elseif ~isempty(trial(i_trial).leftRelTimeTrial) && isempty(trial(i_trial).rightRspTimeTrial)
                        trial(i_trial).SwitchArrival = [];
                        trial(i_trial).SwitchDepart = [];
                    % Use this for the standard switch strategy - silence if using a different strategy
                    elseif isempty(trial(i_trial).leftRelTimeTrial) && ~isempty(trial(i_trial).rightRspTimeTrial)
                        trial(i_trial).SwitchArrival = [];
                        trial(i_trial).SwitchDepart = [];
                    % Use this for the no short response strategy - cautiously interpret this data
    %                 elseif isempty(trial(j).leftRelTimeTrial) && ~isempty(trial(j).rightRspTimeTrial)
    %                     trial(j).SwitchArrival = [];
    %                     trial(j).SwitchDepart = min(trial(j).rightRspTimeTrial(trial(j).rightRspTimeTrial < trial(j).programmedDuration/1000));
                    elseif isempty(trial(i_trial).leftRelTimeTrial) && any(trial(i_trial).rightRspTimeTrial > trial(i_trial).programmedDuration/1000)
                        trial(i_trial).SwitchArrival = [];
                        trial(i_trial).SwitchDepart = [];
                    else
                        trial(i_trial).SwitchArrival = [];
                        trial(i_trial).SwitchDepart = [];
                    end
                    trial(i_trial).ShortRsp = trial(i_trial).leftRspTimeTrial;
                    trial(i_trial).LongRsp = trial(i_trial).rightRspTimeTrial;

                % 1. Start at short and then switch / 2. Start at short and then switch, back to short / 3. Multiple switches
                % 4. Start at long, move to short and then switch / 5. Start and stay at long / 6. Start and stay at short  
                elseif type == 1    % 1 indicates that a short latency trial is rewarded at the right nose poke             
                    if ~isempty(trial(i_trial).rightRelTimeTrial) && ~isempty(trial(i_trial).leftRspTimeTrial) 
                        trial(i_trial).SwitchArrival = min(trial(i_trial).leftRspTimeTrial(trial(i_trial).leftRspTimeTrial > min(trial(i_trial).rightRelTimeTrial)));
                            if ~isempty(trial(i_trial).SwitchArrival) 
                                trial(i_trial).SwitchDepart = max(trial(i_trial).rightRelTimeTrial(trial(i_trial).rightRelTimeTrial < trial(i_trial).SwitchArrival));
                            else trial(i_trial).SwitchArrival=[]; 
                            end
                    elseif isempty(trial(i_trial).rightRelTimeTrial) && isempty(trial(i_trial).leftRspTimeTrial)
                        trial(i_trial).SwitchArrival = [];
                        trial(i_trial).SwitchDepart = [];
                    elseif ~isempty(trial(i_trial).rightRelTimeTrial) && isempty(trial(i_trial).leftRspTimeTrial)
                        trial(i_trial).SwitchArrival = [];
                        trial(i_trial).SwitchDepart = [];
                    % Use this for the standard switch strategy - silence if using a different strategy
                    elseif isempty(trial(i_trial).rightRelTimeTrial) && ~isempty(trial(i_trial).leftRspTimeTrial)
                        trial(i_trial).SwitchArrival = [];
                        trial(i_trial).SwitchDepart = [];
                    % Use this for the no short response strategy - cautiously interpret this data
    %                 elseif isempty(trial(j).rightRelTimeTrial) && ~isempty(trial(j).leftRspTimeTrial)
    %                     trial(j).SwitchArrival = [];
    %                     trial(j).SwitchDepart = min(trial(j).leftRspTimeTrial(trial(j).leftRspTimeTrial < trial(j).programmedDuration/1000);
                    elseif isempty(trial(i_trial).rightRelTimeTrial) && any(trial(i_trial).leftRspTimeTrial > trial(i_trial).programmedDuration/1000)
                        trial(i_trial).SwitchArrival = [];
                        trial(i_trial).SwitchDepart = [];
                    else
                        trial(i_trial).SwitchArrival = [];
                        trial(i_trial).SwitchDepart = [];
                    end
                    trial(i_trial).ShortRsp = trial(i_trial).rightRspTimeTrial;
                    trial(i_trial).LongRsp = trial(i_trial).leftRspTimeTrial;
                else
                    error('')
                end               

                %% response times relative to trial start - **might be able to take this out**
                trial(i_trial).leftResponseTimes = rpInLeft - curTS;
                trial(i_trial).rightResponseTimes = rpInRight - curTS;
                trial(i_trial).backResponseTimes = rpInBack - curTS;
                % release times relative to trial start
                trial(i_trial).leftReleaseTimes = rpOutLeft - curTS;
                trial(i_trial).rightReleaseTimes = rpOutRight - curTS;
                trial(i_trial).backReleaseTimes = rpOutBack - curTS;

                %% time from ITI end to trial initiation poke
                x = rpInBack - trialStart(i_trial);
                trial(i_trial).initiationRT = min(x(x>=0)); % timestamp of first response after the trial was initiated

                %% reaction time: time from signal offset to first response
                startPoke = min(rpInBack(rpInBack >= trialStart(i_trial))); % timestamp of first back poke after trial start
                durationElapsed = startPoke+(trialType(i_trial)./1000); % time that the duration elapsed
                potentialRTs = [min(rpInLeft(rpInLeft>=durationElapsed)) min(rpInRight(rpInRight>=durationElapsed))]; % finds first responses after duration elapsed on left and right nosepokes
                trial(i_trial).RT = min(potentialRTs) - durationElapsed; 

                %% reward entries/exits relative to trial end
                trial(i_trial).rewardEntryTimes = fzIn - curTE;
                trial(i_trial).rewardExitTimes = fzIn - curTE;
                trial(i_trial).reward = reward;
                trial(i_trial).reward_inTrial = reward(reward>curTS & reward<=curTE+0.2) - curTS;

            end
            trial(1).mpc = mpcParsed(lineIdx); % just saves everything
            TrialAnSt = trial; %
        end
    else
    end

end