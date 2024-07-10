function data = readOEphys(dataDir, dataClass, animal_id, oephys_name, dataRange, dataType)
%{
% data = readEphys(dataDir)
% Load OpenEphys Binary and Kwik type continuous, event, or spike data files into Matlab.
% Inputs:
% dataDir: path to data directory
% Outputs:
% data: struct of data with varible field depending on file type and requested data type
%       No data type is converted, it is original data type to avoid redundant conversion and save memory
%       continuous data and spike waveforms need to be converted to be mV unit
%       
% Usage example:
%     dataDir = 'C:\OpenEphys\GUI\2016-08-17_10-18-31';
%     % return data with defalut return data type
%     data = readOEphys(dataDir, 'info'); 
%     data = readOEphys(dataDir, 'events'); 
%     data = readOEphys(dataDir, 'continuous', 1); 
%     data = readOEphys(dataDir, 'spikes', 1); 
% 
%     % get only specfic range for continuous data, useful for very long recording
%     data = readOEphys(dataDir, 'continuous', [1 16; 1 10000]); % return data of first 10000 samples from first 16 channels
% 
%     % get non-default data output
%     data = readOEphys(dataDir, 'continuous', 1, {'ts' 'nsamples' 'data' });
% 
%     % to get available block and data types
%     data = readOEphys('100_CH2.continuous', [1 1], {''}); % no actual data read, but will have header, total block size, available data types
%     dataTypes = {data.dblock.Str}; % extract available data type
%     data = readOEphys('100_CH2.continuous', [1 data.totalBlock-10], dataTypes([1, 3, 4])); % Select data type to read
%     

      
Usage example:

    % return data with defalut return data type
    data = readOEphys('100_CH2.continuous'); 

    % get only specfic data blocks, useful for very long recording
    data = readOEphys('100_CH2.continuous', [1 10]); % return data of 1-10 blocks
    data = readOEphys('100_CH2.continuous', [10 Inf]); % return data of 10-end of blocks

    % get non-default data output
    data = readOEphys('100_CH2.continuous', [1 Inf], {'ts' 'nsamples' 'data' });

    % to get available block and data types
    data = readOEphys('100_CH2.continuous', [1 1], {''}); % no actual data read, but will have header, total block size, available data types
    dataTypes = {data.dblock.Str}; % extract available data type
    data = readOEphys('100_CH2.continuous', [1 data.totalBlock-10], dataTypes([1, 3, 4])); % Select data type to read
    
    % raw continuous data to mV
    contidata = double(data.data) .* data.header.bitVolts;

    % raw spike data to mV
    spikedata = (double(data.data)-32768)./ permute(repmat(data.gain/1000,[1 1 40]), [1 3 2]);

% change log
% 11/2/2017 - fix binary continuos file name, matching file name instead of fixed name (100_ or 101_)

Copyright (C) 2016 Open Ephys

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

<http://www.gnu.org/licenses/>.

%}
% dataDir = 'D:\OpenEphys\GUI\2015-04-06_11-46-58';
% dataDir = 'D:\OpenEphys\GUI\2015-04-06_15-55-16';
% dataDir = 'D:\PlexonData\OpenEphys\K6-FI12_NS-D4_2015-11-23_09-18-26';
% dataDir = 'D:\PlexonData\OpenEphys\K8_2015-12-01_14-27-06';
% dataDir = 'D:\OpenEphys\GUI\2015-08-11_16-37-24';
% dataDir = 'D:\OpenEphys\GUI\2015-08-17_10-18-31';
if exist(dataDir,'dir') == 7
    dirFiles = dir(dataDir);
else
    error('data folder do not exist')
end

% get types of files exist
fileTypes = {'xml' 'spikes' 'continuous' 'events' 'kwe' 'kwx' 'kwd' 'kwik'};
for typeIdx = 1:length(fileTypes)
    fext = sprintf('\\.%s',fileTypes{typeIdx});
    fIdx = ~cellfun(@isempty, regexp({dirFiles.name}, fext, 'once'));
    dfiles.(fileTypes{typeIdx}) = {dirFiles(fIdx).name};
end
dirFiles = dirFiles(~ismember({dirFiles.name},{'.','..'}));
dirFiles = dirFiles(~contains({dirFiles.name},{'._'})); 
pt1 = length(dirFiles) == 1 && dirFiles.isdir == 1 && strncmp(dirFiles.name,'Record Node',11);
pt2 = all(~cellfun(@isempty, regexp({dirFiles.name}, '.xml', 'once')) | (~cellfun(@isempty, regexp({dirFiles.name}, 'exp', 'once')) & [dirFiles.isdir]));

if pt1 || pt2
    if pt1
        dirFiles = dir(fullfile(dirFiles.folder,dirFiles.name));
        dirFiles = dirFiles(~ismember({dirFiles.name},{'.','..'}));
        dirFiles = dirFiles(~contains({dirFiles.name},{'._'})); 
    end
    fIdx = ~cellfun(@isempty, regexp({dirFiles.name}, '.xml', 'once'));
    dfiles.xml = fullfile(dirFiles(fIdx).folder,dirFiles(fIdx).name);
    fIdx = ~cellfun(@isempty, regexp({dirFiles.name}, 'exp', 'once'));
    dirFiles = dir(fullfile(dirFiles(fIdx).folder,dirFiles(fIdx).name));
    dirFiles = dirFiles(~ismember({dirFiles.name},{'.','..'}));
    dirFiles = dirFiles(~contains({dirFiles.name},{'._'}));
    fIdx = ~cellfun(@isempty, regexp({dirFiles.name}, 'rec', 'once'));
    dirFiles = dir(fullfile(dirFiles(fIdx).folder,dirFiles(fIdx).name));    
    fIdx = ~cellfun(@isempty, regexp({dirFiles.name}, '.oebin', 'once'));
    dfiles.oebin = fullfile(dirFiles(fIdx).folder,dirFiles(fIdx).name);
    fIdx = ~cellfun(@isempty, regexp({dirFiles.name}, '^sync_messages.txt', 'once'));
    dfiles.synctxt = fullfile(dirFiles(fIdx).folder,dirFiles(fIdx).name);
end
% decide bin or kwik
if isfield(dfiles,'oebin')
    dataEngine = 'rawbin';  
elseif any([~isempty(dfiles.spikes) ~isempty(dfiles.continuous) ~isempty(dfiles.events)])
    dataEngine = 'bin';
elseif any([~isempty(dfiles.kwe) ~isempty(dfiles.kwx) ~isempty(dfiles.kwd) ~isempty(dfiles.kwik)])
    dataEngine = 'kwik';
else
    error('No data exist')
end

% handling multiple recording in a folder
if strcmp(dataEngine, 'kwik') 
    if numel(dfiles.kwd)>1
        % pick largest kwd
        fsize = zeros(numel(dfiles.kwd),1);
        for i_file = 1:numel(dfiles.kwd)
            fsize(i_file) = dirFiles(strcmp(dfiles.kwd(i_file), {dirFiles.name})).bytes;
        end
        [~,max_idx] = max(fsize);
        % expname = regexp(dfiles.kwd{max_idx},'\w+(?=_\w+)','match','once');
        for typeIdx = 1:length(fileTypes)
            if numel(dfiles.(fileTypes{typeIdx})) == numel(dfiles.kwd)
                dfiles.(fileTypes{typeIdx}) = dfiles.(fileTypes{typeIdx})(max_idx);
            end
        end
    end
end

% decide electrode type and number
if strcmp(dataEngine, 'bin') 
    if ~isempty(dfiles.spikes) && false
        if any(~cellfun(@isempty, regexp(dfiles.spikes, 'SE\w+', 'once')))
            spikeET = 'SE';
        elseif any(~cellfun(@isempty, regexp(dfiles.spikes, 'ST\d+', 'once')))
            spikeET = 'ST';
        elseif any(~cellfun(@isempty, regexp(dfiles.spikes, 'TT\d+', 'once')))
            spikeET = 'TT';
        end
        spkChExist = cellfun(@str2double, regexp(dfiles.spikes, '(?:SE|ST|TT)(\w+).', 'tokens', 'once'));
    end
end

% input check
if ~any(strcmp(dataClass, {'info', 'events', 'continuous', 'spikes'}))
    error('File extension not recognized. Please use a ''.continuous'', ''.spikes'', or ''.events'' file.');
end
if nargin < 6
    dataType = [];
end    
% return data
switch dataClass
    case 'info'
        if ~isempty(dfiles.xml)
            if strcmp(dataEngine, 'rawbin') 
                data = readSessionXML(dfiles.xml);
            else
                data = readSessionXML(fullfile(dataDir,dfiles.xml{1}));
            end
        else
            error('No settings.xml exist')
        end
    case 'events'
        if (nargin > 4) && ~isempty(dataRange), warning('event data do not support channel range and data range'); end
        if strcmp(dataEngine, 'bin') && any(strcmp(dfiles.events,'all_channels.events'))
            data = readOEphysBin(fullfile(dataDir, 'all_channels.events'));
        elseif strcmp(dataEngine, 'kwik') && ~isempty(dfiles.kwe)
            data = readOEphysKwik(fullfile(dataDir, char(dfiles.kwe)));
        elseif strcmp(dataEngine, 'kwik') && ~isempty(dfiles.kwik)
            data = readOEphysKwik(fullfile(dataDir, char(dfiles.kwik)));
        elseif strcmp(dataEngine, 'rawbin') 
            data = readOEphysRawBin(dfiles.oebin,'events',1);
            fileID = fopen(dfiles.synctxt); % extract start time
            line_txt = fgetl(fileID);
            line_txt = fgetl(fileID);
            fclose(fileID); 
            if line_txt == -1 % file is empty for some reason
                %filepath = dataDir;
                cntData = readOEphys(dataDir, 'continuous', animal_id, oephys_name, [1 1; 1 10]);
                %cntData = readOEphys('/Volumes/BovaData1/dLight_data/ASB_rewardmag/ASB4/ASB4_2022-03-08_14-35-55_BASELINE1', 'continuous', [1 1; 1 10]);
                data.startTime = double(cntData.Timestamps(1));
            else 
                data.startTime = str2double(regexp(line_txt,'(?<=start time: )\d+','match','once'));
            end
        else
            error('No event data')
        end
    case 'continuous'
        if nargin > 4
            if numel(dataRange) == 1
                chRange = [dataRange dataRange];
                sampleRange = [1 Inf];
            elseif numel(dataRange) == 4
                chRange = dataRange(1,:);
                sampleRange = dataRange(2,:);
            else
                error('continous data require channel number or channel range/sample range');
            end
        else
            error('continous data require channel number or range');
        end

        if strcmp(dataEngine, 'bin') && ~isempty(dfiles.continuous)
            if diff(chRange) > 0 % multichannel continouse of binary type need to be loop
                chRangeIdx = chRange(1):chRange(2);
                matchedfile = regexp(dfiles.continuous,sprintf('\\d*_CH%d.continuous', chRangeIdx(1)),'match');
                channelfile = matchedfile{~cellfun(@(x) isempty(x),matchedfile)};
                continouseFile = channelfile{1};
                data = readOEphysBin(fullfile(dataDir, continouseFile),sampleRange);
                tempData = zeros(size(data.data,1),length(chRangeIdx),'int16');
                tempData(:,1) = data.data;
                for chIdx = 2:length(chRangeIdx)
                    matchedfile = regexp(dfiles.continuous,sprintf('\\d*_CH%d.continuous', chRangeIdx(chIdx)),'match');
                    channelfile = matchedfile{~cellfun(@(x) isempty(x),matchedfile)};
                    continouseFile = channelfile{1};
                    data = readOEphysBin(fullfile(dataDir, continouseFile),sampleRange);
                    tempData(:,chIdx) = data.data;
                end
                data.data = tempData;
                data.chRange = chRangeIdx;
            else
                matchedfile = regexp(dfiles.continuous,sprintf('\\d*_CH%d.continuous', chRange(1)),'match');
                channelfile = matchedfile{~cellfun(@(x) isempty(x),matchedfile)};
                data = readOEphysBin(fullfile(dataDir, channelfile{1}));
            end
        elseif strcmp(dataEngine, 'kwik') && ~isempty(dfiles.kwd)
            % data = readOEphysKwik(fullfile(dataDir, char(dfiles.kwd)), [chRange; sampleRange]);
            data = readOEphysKwik(fullfile(dataDir, char(dfiles.kwd)), [chRange; sampleRange], dataType);
        elseif strcmp(dataEngine, 'rawbin') 
            data = readOEphysRawBin(dfiles.oebin,'continuous',1,[], [chRange; sampleRange]);
        else
            warning('No continuous data')
            data = [];
        end
    case 'spikes'
        if nargin > 4
            if ~(numel(dataRange) == 1)
                error('spike data do not support channel range, single channel only');
            end
        else
            error('spike data require channel number');
        end
        if strcmp(dataEngine, 'bin') && ~isempty(dfiles.spikes)
            spkFile = sprintf('%s%d.spikes', spikeET, dataRange-1);
            data = readOEphysBin(fullfile(dataDir, spkFile),[1 Inf], dataType);
        elseif strcmp(dataEngine, 'kwik') && ~isempty(dfiles.kwx)
            data = readOEphysKwik(fullfile(dataDir, char(dfiles.kwx)), dataRange, dataType);
        elseif strcmp(dataEngine, 'rawbin') 
            data = load_open_ephys_binary(dfiles.oebin,'spikes',1);            
        else
            error('No spikes data')
        end
end

end


