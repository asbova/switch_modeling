% scriptProcessSingleSortedEphys 
% Processes sorted pl2 files from open ephys recordings to extract waveform 
%
% Data Input:
%
% Output:
%

%%
addpath(genpath('./switch_modeling'))
cd './switch_modeling'

% Identify data directories.
behaviorDataFolder = './data/medpc'; % medPC files
ephysDataFolder = '/Volumes/BovaData1/RawData/modeling_paper'; % open ephys and sorted pl2 files

% Session information - medPC protocols, mouse IDs.
mpcProtocols = {'Switch_18L6R_SITI_RI_MAW', 'Switch_6L18R_SITI_RI_MAW', 'Switch_18L6R_SITI_REINFORCE_FP_V3',...
    'Switch_6L18R_SITI_REINFORCE_FP_V3'}; % List of medpc protocols that were used in these sessions.
ephysMouseFolders = dir(ephysDataFolder); % Each mouse has a sub-folder within the main data folder.
ephysMouseFolders = ephysMouseFolders(~ismember({ephysMouseFolders.name},{'.','..','.DS_Store','._'}));
mouseIDs = {ephysMouseFolders.name};
mpcParsed = getDataIntr(behaviorDataFolder, mpcProtocols, mouseIDs);

%% Create the data structure from directory filenames.

ephysDataStructure = struct;

iCount = 1;
for iMouse = 1 : size(ephysMouseFolders, 1)

    % Find the sorted pl2 files.
    currentMouseFolder = fullfile(ephysDataFolder, ephysMouseFolders(iMouse).name);
    pl2FileNames = dir(currentMouseFolder);
    pl2FileNames = pl2FileNames(contains({pl2FileNames.name},'roughsort') & ~contains({pl2FileNames.name},'._'));
    if isempty(pl2FileNames)
        continue; % No sessions for this mouse have been sorted yet.
    end

    % Identify the open ephys folder(s).
    ephysFileList = dir(currentMouseFolder);
    ephysFileList = ephysFileList(~ismember({ephysFileList.name},{'.', '..'}));
    ephysFileList = ephysFileList([ephysFileList.isdir]);

    for jFile = 1 : length(pl2FileNames)

        % extract mouse ID and session date
        currentFilePathway = fullfile(pl2FileNames(jFile).folder, pl2FileNames(jFile).name);
        [~, currentPL2Name, ~] = fileparts(currentFilePathway);
        currentMouseID = regexp(currentPL2Name, '\w{2,4}\d{1}', 'match', 'once');
        currentDate = regexp(currentPL2Name, '\d+\-\d+\-\d+', 'match', 'once');
        currentDateMPC = datestr(datetime(currentDate, 'InputFormat', 'yyyy-MM-dd'), 'mm/dd/yy');

        ephysDataStructure(iCount).mouseID = currentMouseID;
        ephysDataStructure(iCount).sessionDate = currentDate;
        ephysDataStructure(iCount).mpcDate = currentDateMPC;
        ephysDataStructure(iCount).pl2FilePathway = currentFilePathway;     
        ephysDataStructure(iCount).pl2Name = currentPL2Name;
        ephysDataStructure(iCount).oephysName = ephysFileList(jFile).name;
        ephysDataStructure(iCount).oephysFilePathway = ephysFileList(jFile).folder;
        ephysDataStructure(iCount).format = 'bin'; % DO WE NEED THIS?? FORMAT FOR WHAT??

        iCount = iCount + 1;
    end
end

%% Add behavioral data from MedPC and Open Ephys.

for iSession = 1 : length(ephysDataStructure)
    % MedPC Behavior
    mpcIndex = find(strcmp(ephysDataStructure(iSession).mouseID, {mpcParsed.Subject}) & strcmp(ephysDataStructure(iSession).mpcDate, {mpcParsed.StartDate}));
    if isscalar(mpcIndex)
        ephysDataStructure(iSession).mpcData = mpcParsed(mpcIndex);
        ephysDataStructure(iSession).mpcTrialData = getTrialDataSwitch(ephysDataStructure(iSession).mpcData);
    else
        % Do not analyze data
    end

    % Open Ephys Behavior
    currentEphysPathway = fullfile(ephysDataStructure(iSession).oephysFilePathway, ephysDataStructure(iSession).oephysName);
    ephysDataStructure(iSession).ephysEvents = get_mpc_bin_event_oe3(currentEphysPathway, ephysDataStructure(iSession).mpcData, ...
        ephysDataStructure(iSession).mouseID, ephysDataStructure(iSession).oephysName);   
    [ephysDataStructure(iSession).withinTrialEphysEvents, ephysDataStructure(iSession).trialSpecificEphysEvents] = extractTrialEventsEphys(ephysDataStructure(iSession));
end

%% Neuron analysis

for iSession = 1 : length(ephysDataStructure)
    ephysDataStructure(iSession).neurons = extractNeuronData(ephysDataStructure(iSession));
end

% Clustering analysis
allNeurons = [];
for iSession = 1 : length(ephysDataStructure)
    allNeurons = [allNeurons ephysDataStructure(iSession).neurons];
end
neurons = clusterStriatalNeurons(allNeurons, 0.85);

% Distribute neuron classification identifiers to each session within ephysDataStructure.
startIndex = 1;
totalUnits = 0;
for iSession = 1 : length(ephysDataStructure)
    nUnits = size(ephysDataStructure(iSession).neurons, 2);
    [ephysDataStructure(iSession).neurons.type] = deal(neurons(startIndex : startIndex + nUnits-1).type);
    totalUnits = totalUnits + nUnits;
    startIndex = totalUnits + 1;
end

% Save file.
save('./data/matfiles/neuronStructure.mat', 'ephysDataStructure');