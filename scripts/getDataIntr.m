% getDataIntr.m, Youngcho Kim, Narayanan's Lab 2014
% Interactively get MedPC data parsed with mpcParser.m
% Parses MedPC files of requested file, retrun a structure with all MPC data
% Usage: getDataIntr(foldername); 
%
% the parseMPC is the return structure with metadata(date, subject, MSN, sbjIDs)
% and timing data (S,H, etc)

% sbjIDs can be cell array with exact animal ID or character array with regexp
% mpcParsed = getDataIntr(local_data, MSN, 'RT\d*');
% mpcParsed = getDataIntr(local_data, MSN, {'RT1' 'RT2'});

% History
% 1/9/14  : quickly parse for MSN only, ask what to analyze.
% 1/19/14 : allow to select multiple MSN
% 1/23/14 : add error messages
% 9/25/17 : add subject input option
% Get data folder, then MPC file list
% find the matching MSN, put them in a structure

function mpcParsed = getDataIntr(folder_name, MSN, sbjIDs, dateRange)

if nargin < 3
    sbjIDs = [];
end
if nargin < 4
    dateRange = [];
end

if nargin == 3 && isstruct(sbjIDs)
    group = sbjIDs;
    groupFields = fieldnames(group);
    numOfGroup = size(groupFields,1);
    dayAnimalIdx = 1; 
    for groupIdx = 1:numOfGroup
        dayNanimalGroup = group.(char(groupFields(groupIdx)));
        for i_dayAnimal = 1:size(dayNanimalGroup,2)
            dayNanimal = dayNanimalGroup{i_dayAnimal};
            groupDataTableIdx{dayAnimalIdx,1} = dayNanimal{1};
            groupDataTableIdx{dayAnimalIdx,2} = dayNanimal{2};
            dayAnimalIdx = dayAnimalIdx + 1;
        end
    end
    dateRange = groupDataTableIdx(:,1);
    sbjIDs = groupDataTableIdx(:,2);
    is_group = true;
else
    is_group = false;
end

if isfile(folder_name)
    [folder_name,fileList{1},~] = fileparts(folder_name);
elseif isfolder(folder_name)
    % get MEDPC file list from current directory by matching !2013-01-13
    tmp_fileList = dir(folder_name);
    fileList = regexpi({tmp_fileList.name},'(^\d{4}-\d{2}-\d{2}.txt|^!\d{4}-\d{2}-\d{2}.txt|^!\d{4}-\d{2}-\d{2})','match');
    fileList = [fileList{:}];

    if size(fileList,2)==0
        error('No Matching MEDPC DATA file found!')
    end
end

if iscell(dateRange) && ~isempty(dateRange) % date range is requested
    % date_string = cellfun(@(x) x(2:end), fileList,'UniformOutput', false);
    date_string = regexpi(fileList,'(\d{4}-\d{2}-\d{2})','match','once');
    if length(dateRange) < 3 % 1 or 2 is considered as range
        file_dates = datetime(date_string,'InputFormat','yyyy-MM-dd');
        dateStart = datetime(dateRange{1},'InputFormat','yyyy-MM-dd');
        if length(dateRange)>1
            dateEnd = datetime(dateRange{2},'InputFormat','yyyy-MM-dd');
        elseif length(dateRange)==1
            dateEnd = max(file_dates);
        else

        end
        dateIdx = file_dates>=dateStart & file_dates<=dateEnd;
    else  % more than 2 will be considered as individual dates
        %dateRange = convertDateStr(dateRange,'file_date');
        dateRange = convertDateStr(dateRange,'data_date');
        X = cellfun(@(c)strcmp(c,date_string),dateRange,'UniformOutput',false);
        dateIdx = any(vertcat(X{:}));
    end
    fileList = fileList(dateIdx);
end

% quickly parse MSN only from file list
msnListNum = 1;
msnListByFile = {};
sbjListByFile = {};
for filenum=1:length(fileList)
    fString = fileread(fullfile(folder_name,char(fileList(filenum)))); % read the file as one string
    msns = regexp(fString,'MSN:\s([\w-]+)\r','tokens');
    subjects = regexp(fString,'Subject:\s([\w-\s]*)\r','tokens');
    msnListByFileNum = 1;
    for msnNum = 1:min(length(msns), length(subjects))
        msnListByFile(filenum,msnListByFileNum) = msns{msnNum};
        msnList(msnListNum) = msns{msnNum};
        sbjListByFile(filenum,msnListByFileNum) = subjects{msnNum};
        msnListNum = msnListNum + 1;
        msnListByFileNum = msnListByFileNum + 1;
    end
end


if size(msnListByFile,1)==0
    error('No MEDPC DATA found from files!')
end

if ~exist('MSN','var') || isempty(MSN)
    % get all experiment protocol (MSN) from data, ask what to analyze
    % sort by frequency before ask 
    [uniqueMSN, ~, j] = unique(msnList);
    % if only one MSN is available, don't ask
    if size(uniqueMSN,2) == 1
        MSN = uniqueMSN(1);
    else
        [~, sorted_locations] = sort(hist(j, 1:max(j)));
        uniqueMSN = uniqueMSN(fliplr(sorted_locations)).';
        [Selection,~] = listdlg('PromptString','Select a Protocol', 'ListString',uniqueMSN);
        MSN = uniqueMSN(Selection);
    end
end

% handle regexp subject IDs
if ~iscell(sbjIDs) && ischar(sbjIDs) && ~isempty(sbjIDs)
    matchSbj = regexp(sbjListByFile, sbjIDs,'match');
    sbjIDs = unique([matchSbj{:}]);
end

if iscell(sbjIDs)&& ~isempty(sbjIDs)
    % find file list that have requested MSN & subject IDs
    for msnSelected = 1:length(MSN)
        for filenum=1:length(fileList)
            tmp_msnbyFile = msnListByFile(filenum,:);
            tmp_msnbyFile = tmp_msnbyFile(~cellfun('isempty',tmp_msnbyFile));
            tmp_sbjbyFile = sbjListByFile(filenum,:);
            tmp_sbjbyFile = tmp_sbjbyFile(~cellfun('isempty',tmp_sbjbyFile));
            sbjIn = false;
            for sbjIdx = 1:length(sbjIDs)
                sbjIn = sbjIn || any(strcmp(sbjIDs(sbjIdx),tmp_sbjbyFile));
            end
            parsingFiles(msnSelected,filenum) = sbjIn && any(strcmp(MSN{msnSelected},tmp_msnbyFile));
        end
    end
else % find file list that have requested MSN
    for msnSelected = 1:length(MSN)
        for filenum=1:length(fileList)
            tmp_msnbyFile = msnListByFile(filenum,:);
            tmp_msnbyFile = tmp_msnbyFile(~cellfun('isempty',tmp_msnbyFile));
            parsingFiles(msnSelected,filenum) = ismember(MSN{msnSelected},tmp_msnbyFile);
        end
    end
end

% parse mpc files for requested MSN
mpcParsed = struct;
pfileList = fileList(any(parsingFiles,1));
for filenum=1:length(pfileList)
    fNAME = fullfile(folder_name,char(pfileList(filenum)));
    mpcParsed = mpcParser(fNAME, MSN, mpcParsed, sbjIDs);
end

% remove seesions is not in group var
if is_group
    is_in_group = false(size(groupDataTableIdx,1),1);
    for i_gr = 1:size(groupDataTableIdx,1)
        is_in_group(strcmp(groupDataTableIdx(i_gr,2),{mpcParsed.Subject}) & strcmp(groupDataTableIdx(i_gr,1),{mpcParsed.StartDate})) = true;
    end
    mpcParsed = mpcParsed(is_in_group); 
end

% % remove fields with no data 
% parsedFields = fieldnames(mpcParsed);
% for fieldNum = 1:length(parsedFields)
%     emptyData = 0;
%     for dataNum = 1:size(mpcParsed,2)
%         emptyData = emptyData + isempty(mpcParsed(dataNum).(parsedFields{fieldNum}));
%     end
%     if emptyData == size(mpcParsed,2)
%         mpcParsed = rmfield(mpcParsed,parsedFields{fieldNum});
%     end
% end
% 
% % data check for repeated testing in a day
% uniqueSbj = unique({mpcParsed.Subject});
% for i = 1:length(uniqueSbj)
%    subIdx = strcmp(uniqueSbj(i),{mpcParsed.Subject});
%    dateList = {mpcParsed(subIdx).StartDate};
%    if size(dateList,2) > size(unique(dateList),2)
%        warning('duplicated Animal ID %s found, check data!', char(uniqueSbj(i)));
%    end
% end

% remove data with no subject name
mpcParsed = mpcParsed(~strcmp('0',{mpcParsed.Subject}));

% fix Animal ID
for i=1:size(mpcParsed,2)
    sbj = mpcParsed(i).Subject;
    % check start with letter
    if ~isstrprop(sbj(1), 'alpha')
        sbj = strcat('FIX',sbj);
    end
    % replace "-" with "_"
    sbj = regexprep(sbj,'-','_');
    mpcParsed(i).Subject = sbj;
end

% Sort structure by MSN, subject name and timestamp(start date)
[~, ind] = sortrows([{mpcParsed.MSN}',{mpcParsed.Subject}',{mpcParsed.timeStamp}'], [1 2 3]);
mpcParsed = mpcParsed(ind);

