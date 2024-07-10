function D = readOEphysRawBin(jsonFile, type, index, varargin)
%function D=load_open_ephys_binary(oebinFile, type, index)
%
%Loads data recorded by Open Ephys in Binary format
%  oebinFile: The path for the structure.oebin json file
%  type: The type of data to load. Can be 'continuous', 'events' or
%  'spikes'
%  index: The index of the recorded element as appears in the oebin file
%(See also list_open_ephys_binary to extract recorded elements and indexes)
%
%Returns a structure with the header and relevant loaded data
%
%Example:
%D=load_open_ephys_binary('recording_1/structure.oebin','spikes',2)
%
%When loading continuous data, an optional fourth argument 'mmap' can be
%used:
%D=load_open_ephys_binary('recording_1/structure.oebin','continuous',1,'mmap')
%In this case, the Data member from the returned structure contains not an
%array with the data itself but a memory-mapped object of the file, which
%can be used to access its contents. This helps loading big-sized files, as
%it does not require loading the entire file in memory.
%In this case, the data may be accessed using the field D.Data.Data(1).mapped
%For example: D.Data.Data(1).mapped(chan,startSample:endSample)
%
%
%Limitations:
%-TEXT events are not supported by the NPY reading library, so their data
%will not load, but all other fields will.
%-Some metadata structures might not be supported by the NPY reading
%library. In this case the metadata will not be loaded
%In both cases a warning message will be displayed but the program will not
%fail
%
%
%This functions requires the functions readNPY and readNPYHeader 
%from npy-matlab package from kwik-team
%(https://github.com/kwikteam/npy-matlab)
%Requires minimum MATLAB version: R2016b
% 
% if (exist('readNPY.m','file') == 0)
%     error('OpenEphys:loadBinary:npyLibraryNotfound','readNPY not found. Please be sure that npy-matlab is accessible');
% end
% 
% if (exist('readNPYheader.m','file') == 0)
%     error('OpenEphys:loadBinary:npyLibraryNotfound','readNPYheader not found. Please be sure that npy-matlab is accessible');
% end
if (nargin > 3 && strcmp(varargin{1},'mmap'))
    continuousmap = true;
else
    continuousmap = false;
end
if (nargin > 4)
    dRange = varargin{2};
else
    dRange = [];
end


json=jsondecode(fileread(jsonFile));

%Load appropriate header data from json
switch type
    case 'continuous'
        header=json.continuous(index);
    case 'spikes'
        header=json.spikes(index);
    case 'events'
        if (iscell(json.events))
            header=json.events{index}; 
        else
            header=json.events(index);
        end
    otherwise
        error('Data type not supported');
end

%Look for folder
f=java.io.File(jsonFile);
if (~f.isAbsolute())
    f=java.io.File(fullfile(pwd,jsonFile));
end
f=java.io.File(f.getParentFile(),fullfile(type, header.folder_name));
if(~f.exists())
    error('Data folder not found');
end
folder = char(f.getCanonicalPath());
D=struct();
D.Header = header;

switch type
    case 'continuous'
        %D.Timestamps = readNPY(fullfile(folder,'timestamps.npy'));
        D.Timestamps = readNPY_TS(fullfile(folder,'timestamps.npy'));
        contFile=fullfile(folder,'continuous.dat');
        if isempty(dRange)
            if (continuousmap)
                file=dir(contFile);
                samples=file.bytes/2/header.num_channels;
                D.Data=memmapfile(contFile,'Format',{'int16' [header.num_channels samples] 'mapped'});
            else
                file=fopen(contFile);
                D.Data=fread(file,[header.num_channels Inf],'int16');
                fclose(file);
            end
        else
            file = dir(contFile);
            sampleBytes = 2;
            totalSampleNum = file.bytes/sampleBytes/header.num_channels;
            if numel(dRange) == 1
                chRange = [dRange(1) dRange(1)];
                sampleRange = [1 totalSampleNum];
            elseif dRange(2,2) == Inf
                chRange = dRange(1,:);
                sampleRange = [1 totalSampleNum];
            elseif numel(dRange) == 4
                chRange = dRange(1,:);
                sampleRange = dRange(2,:);
            end
            chNum = header.num_channels;
            sampleNum = sampleRange(2) - sampleRange(1) +1;
            fid = fopen(contFile);
            fseek(fid, header.num_channels * (sampleRange(1)-1)* sampleBytes, 'bof'); % sampleRange(1) CH1
            D.data = fread(fid, [chNum sampleNum], sprintf('%d*int16=>int16',chNum), header.num_channels * sampleBytes - chNum*sampleBytes, 'l');
            D.data = D.data(chRange(1):chRange(2),:)';
            fclose(fid);
            if sampleRange(2) > size(D.Timestamps,1)
                sampleRange(2) = size(D.Timestamps,1);
            end
            D.Timestamps = D.Timestamps(sampleRange(1):sampleRange(2));
        end
    case 'spikes'
        D.Timestamps = readNPY(fullfile(folder,'spike_times.npy'));
        D.Waveforms = readNPY(fullfile(folder,'spike_waveforms.npy'));
        D.ElectrodeIndexes = readNPY(fullfile(folder,'spike_electrode_indices.npy'));
        D.SortedIndexes = readNPY(fullfile(folder,'spike_clusters.npy'));
    case 'events'
        D.Timestamps = readNPY(fullfile(folder,'timestamps.npy'));
        D.ChannelIndex = readNPY(fullfile(folder,'channels.npy'));
        f=java.io.File(folder);
        group=char(f.getName());
        if (strncmp(group,'TEXT',4))
            %D.Data = readNPY(fullfile(folder,'text.npy'));
            warning('TEXT files not supported by npy library');
        elseif (strncmp(group,'TTL',3))
            D.Data = readNPY(fullfile(folder,'channel_states.npy'));
            wordfile = fullfile(folder,'full_words.npy');
            if (isfile(wordfile))
                D.FullWords = readNPY(wordfile);
            end
        elseif (strncmp(group,'BINARY',6))
           D.Data = readNPY(fullfile(folder,'data_array.npy'));
        end       
end
metadatafile = fullfile(folder,'metadata.npy');
if (isfile(metadatafile))
    try
    D.MetaData = readNPY(metadatafile);
    catch EX
        fprintf('WARNING: cannot read metadata file.\nData structure might not be supported.\n\nError message: %s\nTrace:\n',EX.message);
        for i=1:length(EX.stack)
            fprintf('File: %s Function: %s Line: %d\n',EX.stack(i).file,EX.stack(i).name,EX.stack(i).line);
        end
    end
end
end

function data = readNPY_TS(filename)
    persistent TS_index;
    if isempty(TS_index)
        TS_index = {};
        TS_index{1}.Data = readNPY(filename);
        TS_index{1}.Name = filename;
        data = TS_index{1}.Data;
        return;
    else
        if strcmp(TS_index{1}.Name, filename)
            if isfield(TS_index{1}, 'Data')
                % if the index structure has pl2 index assigned, return it
                data = TS_index{1}.Data;
                return;
            else
                TS_index = {};
                TS_index{1}.Data = readNPY(filename);
                TS_index{1}.Name = filename;
                data = TS_index{1}.Data;
                return;
            end
        else
            TS_index = {};
            TS_index{1}.Data = readNPY(filename);
            TS_index{1}.Name = filename;
            data = TS_index{1}.Data;
            return;
        end

    end    
end

function data = readNPY(filename)
% Function to read NPY files into matlab.
% *** Only reads a subset of all possible NPY files, specifically N-D arrays of certain data types.
% See https://github.com/kwikteam/npy-matlab/blob/master/tests/npy.ipynb for
% more.
%

[shape, dataType, fortranOrder, littleEndian, totalHeaderLength, ~] = readNPYheader(filename);

if littleEndian
    fid = fopen(filename, 'r', 'l');
else
    fid = fopen(filename, 'r', 'b');
end

try

    [~] = fread(fid, totalHeaderLength, 'uint8');

    % read the data
    data = fread(fid, prod(shape), [dataType '=>' dataType]);

    if length(shape)>1 && ~fortranOrder
        data = reshape(data, shape(end:-1:1));
        data = permute(data, [length(shape):-1:1]);
    elseif length(shape)>1
        data = reshape(data, shape);
    end

    fclose(fid);

catch me
    fclose(fid);
    rethrow(me);
end
end


function [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename)
% function [arrayShape, dataType, fortranOrder, littleEndian, ...
%       totalHeaderLength, npyVersion] = readNPYheader(filename)
%
% parse the header of a .npy file and return all the info contained
% therein.
%
% Based on spec at http://docs.scipy.org/doc/numpy-dev/neps/npy-format.html

fid = fopen(filename);

% verify that the file exists
if (fid == -1)
    if ~isempty(dir(filename))
        error('Permission denied: %s', filename);
    else
        error('File not found: %s', filename);
    end
end

try
    
    dtypesMatlab = {'uint8','uint16','uint32','uint64','int8','int16','int32','int64','single','double', 'logical'};
    dtypesNPY = {'u1', 'u2', 'u4', 'u8', 'i1', 'i2', 'i4', 'i8', 'f4', 'f8', 'b1'};
    
    
    magicString = fread(fid, [1 6], 'uint8=>uint8');
    
    if ~all(magicString == [147,78,85,77,80,89])
        error('readNPY:NotNUMPYFile', 'Error: This file does not appear to be NUMPY format based on the header.');
    end
    
    majorVersion = fread(fid, [1 1], 'uint8=>uint8');
    minorVersion = fread(fid, [1 1], 'uint8=>uint8');
    
    npyVersion = [majorVersion minorVersion];
    
    headerLength = fread(fid, [1 1], 'uint16=>uint16');
    
    totalHeaderLength = 10+headerLength;
    
    arrayFormat = fread(fid, [1 headerLength], 'char=>char');
    
    % to interpret the array format info, we make some fairly strict
    % assumptions about its format...
    
    r = regexp(arrayFormat, '''descr''\s*:\s*''(.*?)''', 'tokens');
    dtNPY = r{1}{1};    
    
    littleEndian = ~strcmp(dtNPY(1), '>');
    
    dataType = dtypesMatlab{strcmp(dtNPY(2:3), dtypesNPY)};
        
    r = regexp(arrayFormat, '''fortran_order''\s*:\s*(\w+)', 'tokens');
    fortranOrder = strcmp(r{1}{1}, 'True');
    
    r = regexp(arrayFormat, '''shape''\s*:\s*\((.*?)\)', 'tokens');
    shapeStr = r{1}{1}; 
    arrayShape = str2num(shapeStr(shapeStr~='L'));

    
    fclose(fid);
    
catch me
    fclose(fid);
    rethrow(me);
end
end