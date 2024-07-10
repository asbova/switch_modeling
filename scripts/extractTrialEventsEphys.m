function [trialEventStructure, eventStructure] = extractTrialEventsEphys(sessionData)

% Pulls out within trial events for each session
%
% Input:
%   sessionData:    Structure containing open ephys and medpc beahavioral data for the current session
%
% Output:
%   

    % Session information
    if contains(sessionData.mpcData.MSN, '6L18R')
        protocol = 0; % 6L18R
    else
        protocol = 1; % 18L6R
    end
    trialDurations = [sessionData.mpcTrialData.programmedDuration];
    trialOutcome = [sessionData.mpcTrialData.outcome]; % 1 = incorrect correction trial; 2 = correct correction trial;
                                                       % 3 = incorrect trial; 4 = correct trial

    % Classify trials based on trial duration and outcome, from MedPC.
    trialsWithSwitch = ~cellfun(@isempty, {sessionData.mpcTrialData.SwitchDepart});
    longTrialsWithSwitch = (trialsWithSwitch == 1) & (trialDurations == 18000);
    correctShortTrials = (trialDurations == 6000 & (trialOutcome == 4 | trialOutcome == 2));
    correctLongTrials = (trialDurations == 18000 & (trialOutcome == 4 | trialOutcome == 2));
    incorrectShortTrials = (trialDurations == 6000 & (trialOutcome == 3 | trialOutcome == 1));
    incorrectLongTrials = (trialDurations == 18000 & (trialOutcome == 3 | trialOutcome == 1));
    
    % Extract Open Ephys events.
    cueOnTimestamps = sessionData.ephysEvents.evt3.ts;
    trialEndTimestamps = sessionData.ephysEvents.evt17.ts; % Not cues off.
    leftResponseTimestamps = sessionData.ephysEvents.evt7.ts;
    leftReleaseTimestamps = sessionData.ephysEvents.evt11.ts;
    rightResponseTimestamps = sessionData.ephysEvents.evt19.ts;
    rightReleaseTimestamps = sessionData.ephysEvents.evt15.ts;
    rewardDispenseTimestamps = sessionData.ephysEvents.evt9.ts;

    % % Remove incorrect trials. Doesn't seem necessary...?
    % if ~isempty(find(diff(cueOnTimestamps)<0.001))
    %     badTimestamps = find(diff(cueOnTimestamps)<0.001) + 1;
    %     correctShortTrials(badTimestamps) = [];
    %     correctLongTrials(badTimestamps) = [];
    %     incorrectShortTrials(badTimestamps) = [];
    %     incorrectLongTrials(badTimestamps) = [];
    %     longTrialsWithSwitch(badTimestamps) = [];
    %     trialEndTimestamps(badTimestamps) = [];
    %     cueOnTimestamps(badTimestamps) = [];
    % else
    %     % 
    % end

    % NEED TO ADD CODE TO FIGURE OUT IF DIFFERENCE BETWEEN CUE ON AND TRIAL END TIMESTAMPS IS ERROR?
    nTrials = min(length(cueOnTimestamps), length(trialEndTimestamps));
    trialEventStructure = struct;

    for iTrial = 1 : nTrials

        currentTrialStart = cueOnTimestamps(iTrial);
        currentTrialEnd = trialEndTimestamps(iTrial);
        trialEventStructure(iTrial).cuesOnTime = currentTrialStart;

        % Extract timestamps of nosepoke events that are within the trial.
        trialEventStructure(iTrial).leftResponseTime = leftResponseTimestamps(leftResponseTimestamps > currentTrialStart & leftResponseTimestamps <= currentTrialEnd) - currentTrialStart;
        trialEventStructure(iTrial).leftReleaseTime = leftReleaseTimestamps(leftReleaseTimestamps > currentTrialStart & leftReleaseTimestamps <= currentTrialEnd) - currentTrialStart;
        trialEventStructure(iTrial).rightResponseTime = rightResponseTimestamps(rightResponseTimestamps > currentTrialStart & rightResponseTimestamps <= currentTrialEnd) - currentTrialStart;
        trialEventStructure(iTrial).rightReleaseTime = rightReleaseTimestamps(rightReleaseTimestamps > currentTrialStart & rightReleaseTimestamps <= currentTrialEnd) - currentTrialStart;

        % Identify timestamps of switch responses.
        if protocol == 0 && longTrialsWithSwitch(iTrial) == 1
            validLongResponses = trialEventStructure(iTrial).rightResponseTime > min(trialEventStructure(iTrial).leftReleaseTime);
            trialEventStructure(iTrial).switchArrivalTime = min(trialEventStructure(iTrial).rightResponseTime(validLongResponses));
            validShortRespones = trialEventStructure(iTrial).leftReleaseTime < trialEventStructure(iTrial).switchArrivalTime;
            trialEventStructure(iTrial).switchDepartureTime = max(trialEventStructure(iTrial).leftReleaseTime(validShortRespones));
        elseif protocol == 1 && longTrialsWithSwitch(iTrial) == 1
            validLongResponses = trialEventStructure(iTrial).leftResponseTime > min(trialEventStructure(iTrial).rightReleaseTime);
            trialEventStructure(iTrial).switchArrivalTime = min(trialEventStructure(iTrial).leftResponseTime(validLongResponses));
            validShortRespones = trialEventStructure(iTrial).rightReleaseTime < trialEventStructure(iTrial).switchArrivalTime;
            trialEventStructure(iTrial).switchDepartureTime = max(trialEventStructure(iTrial).rightReleaseTime(validShortRespones));
        else
            % There is no switch and/or this is not a long trial.
        end
    end

    eventStructure = struct;
    eventStructure.correctShortTrialCuesOn = cueOnTimestamps(correctShortTrials);
    eventStructure.correctLongTrialCuesOn = cueOnTimestamps(correctLongTrials);
    eventStructure.incorrectShortTrialCuesOn = cueOnTimestamps(incorrectShortTrials);
    eventStructure.incorrectLongTrialCuesOn = cueOnTimestamps(incorrectLongTrials);
    eventStructure.trialsWithSwitchCuesOn = cueOnTimestamps(longTrialsWithSwitch);
    eventStructure.switchDeparture = [trialEventStructure(longTrialsWithSwitch).cuesOnTime] + [trialEventStructure(longTrialsWithSwitch).switchDepartureTime];
    eventStructure.switchArrival = [trialEventStructure(longTrialsWithSwitch).cuesOnTime] + [trialEventStructure(longTrialsWithSwitch).switchArrivalTime];
    eventStructure.rewardDispense = rewardDispenseTimestamps;
    if protocol == 0
        eventStructure.shortResponse = leftResponseTimestamps;
        eventStructure.longResponse = rightResponseTimestamps;
    elseif protocol == 1
        eventStructure.shortResponse = rightResponseTimestamps;
        eventStructure.longResponse = leftResponseTimestamps;
    end


end