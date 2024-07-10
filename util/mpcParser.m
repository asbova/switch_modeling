% mpcParser.m, Youngcho Kim, Narayanan's Lab 2014
% Parses MedPC files of requested file, retrun a structure with matched MPC data
% Usage: mpcParser(filename, MSN, output-structure); 
%
% the parseMPC is the return structure with metadata(date, subject, MSN)
% and timing data (S,H, etc)

% History
% 12/7/13 : rewrite most parts to automatically parse whole header and
% data. 5x speed improvement by removing 0.000 before str2double
% 1/9/14  : add input for MSN to parse only selected MSN's data. improved
% efficiency of parsing
% 1/23/14 : change Segmentation of experiments code for accurate MSN match
% 7/18/14 : accurately remove trailing zeros from data
% 1/23/15 : change time stamp to include time (h,m,s)
% 1/27/15 : take cellarray MSN
% 4/20/15 : use sscanf for str2double (3-4 times faster)
% 12/18/15 : merge into big loop, remove a redundant regex. 2 times faster
% 9/25/17 : add subject input option
% following are data format from MEDPC output file

%{
Start Date: 11/15/13
End Date: 11/15/13
Subject: C8M1
Experiment: NOTONE D1
Group: 0
Box: 6
Start Time:  5:38:54
End Time:  6:34:32
MSN: MOUSE_FI12_notone_MS
A:      10.000
G:
     0:       -0.025       -0.020       -0.015       -0.010       -0.005
     5:        0.005        0.010        0.015        0.020        0.025
%}

%{ 
debug
parse a MEDPC file, return structure
clear all

fName = '!2014-05-07';
mpcParsed = struct;
MSN = 'MOUSE_FI12';

fNAME =fullfile(folder_name,fNAME)
%}

% parse a MEDPC file, return structure
function mpcParsed = mpcParser(fName, MSN, mpcParsed, sbjIDs)

if nargin < 4
    sbjIDs = [];
end

% Segment each experiment into cells with a keyword (Start Date:)
% store them, if matching MSN exist
fString = fileread(fName); 
expsStartPts = strfind(fString, 'Start Date:'); % Case sensitive
msns = regexp(fString,'MSN:\s([\w-]+)\r','tokens');
matchMSNIdx = cellfun(@any, cellfun(@(x) strcmp(x,MSN), [msns{:}],'UniformOutput',false));
if sum(matchMSNIdx)==0; error('No Data with the MSN!'); end

if iscell(sbjIDs)&& ~isempty(sbjIDs)
    subjects = regexp(fString,'Subject:\s([\w-\s]*)\r','tokens');
    matchSbjIdx = cellfun(@any, cellfun(@(x) strcmp(x,sbjIDs), [subjects{:}],'UniformOutput',false));
else
    matchSbjIdx = true(size(matchMSNIdx));
end
if length(matchMSNIdx) ~= length(matchSbjIdx)
    return
end
matchMSNIdx = matchMSNIdx & matchSbjIdx;

startPts = expsStartPts(matchMSNIdx); 
endPtIdx = find(matchMSNIdx)+1;
endPts = expsStartPts(endPtIdx(endPtIdx <= length(expsStartPts)))-1; 
if sum(endPtIdx>length(expsStartPts))==1; endPts(end+1)= size(fString,2); end
if isempty(fieldnames(mpcParsed)), ctsize = 0; else ctsize = size(mpcParsed,2); end
formatIn = 'mm/dd/yyyy HH:MM:SS';
% Segment Header and Data parts
% Segments Data into Each sub-Data (A to Z)
for i = 1:length(startPts)
    headerEnd = regexp(fString(startPts(i):endPts(i)),'Start Date:(.+)MSN:\s[\w-]+\r','end');
    HeaderString = fString(startPts(i):startPts(i)+headerEnd-1);
    headerLines = regexp(HeaderString,'\n','split');
    for j = 1:length(headerLines)
        header = regexp(headerLines(j),':','split','once');
        mpcParsed(ctsize+i).(sscanf(header{1}{1},'%s')) = strtrim(header{1}{2});
    end
    mpcParsed(ctsize+i).timeStamp = datenum([mpcParsed(ctsize+i).StartDate ' ' mpcParsed(ctsize+i).StartTime],formatIn);
    dataString = fString(startPts(i)+headerEnd:endPts(i));
    dataStartPt = regexp(dataString,'[A-Z]:','start');
    dataEndPt = [dataStartPt(2:end)-1 length(dataString)];
    for j = 1:length(dataStartPt)
        subData = dataString(dataStartPt(j):dataEndPt(j));
        data_match = regexp(subData,'\d{1,10}\.\d{3}','match'); % literal match on .
        dataDouble = sscanf(sprintf('%s,', data_match{:}), '%f,'); % 3-4 times faster than str2double
        mpcParsed(ctsize+i).(regexp(subData,'[A-Z]','match','once')) = dataDouble(1:find(dataDouble,1,'last'));
    end
end

