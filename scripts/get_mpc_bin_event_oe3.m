function events = get_mpc_bin_event_oe3(event_datapath, mpcParsed, animal_id, oephys_name)

% Binary event
event_record_type = {...
'1 0 0 0 0', 'SESSION START'; ...
'0 1 0 0 0', 'Back Nosepoke Light ON'; ...
'0 0 1 0 0', 'Back Nosepoke Light OFF'; ...
'0 0 0 1 0', 'ITI START'; ...
'0 0 0 0 1', 'ITI END'; ...
'1 1 0 0 0', 'CUES ON'; ...
'1 0 1 0 0', 'CUES OFF'; ...
'1 0 0 1 0', 'REWARD DISPENSE'; ...
'1 0 0 0 1', 'TRIAL END'; ...
'1 1 1 0 0', 'LEFT RESPONSE'; ...
'1 1 0 1 0', 'LEFT RELEASE'; ...
'1 1 0 0 1', 'RIGHT RESPONSE'; ...
'1 1 1 1 0', 'RIGHT RELEASE'; ...
'1 1 1 1 1', 'BACK RESPONSE'; ...
'0 1 1 0 0', 'BACK RELEASE'; ...
'0 1 0 1 0', 'REWARD RESPONSE/CHECKING FEEDER'; ...
'0 1 0 0 1', 'REWARD RELEASE/REMOVE HEAD FROM FEEDER';... 
'0 0 1 0 1', 'SHORT TRIAL';...
'1 0 1 0 1', 'LONG TRIAL';...
'0 1 1 0 1', 'NO STIM TRIAL';...
'1 1 1 0 1', 'STIM TRIAL'...
};

% Z^zRPBreakBACK; Z^zPlexRESPBACK; -> Z^zPlexCUESON;, Z^zPlexNPLBACKOFF;

% convert binary event to decimal
event_type_dec = cellfun(@(x) bin2dec(fliplr(x)), event_record_type(:,1));
event_record_type(1:21,3) = num2cell(event_type_dec);

% MPC -> Superport -> OE 
% C5  -> 16 -> 0 
% TTL -> 15 -> 1
% C4  -> 14 -> 2 
% C3  -> 13 -> 3 
% C2  -> 12 -> 4
% C1  -> 11 -> 5 
eventData = readOEphys(event_datapath, 'events', animal_id, oephys_name);
[~,name,~] = fileparts(event_datapath);
fprintf('\n')
fprintf('%s \n',name)
fprintf('\n')

% 11/1/2021: fix event start time from raw open-ephys format
startTick = eventData.startTime; % re-align event time stamps to be 0 at processor 100 (FPGA buffer time)
evt_ts = double(eventData.Timestamps) - startTick;
evt_data = double(eventData.Data);

% check glitch, single channel rise/fall in 1 tick
evt_debug = [];
for ch_num = [1 3:6]
    fall_idx = find(eventData.Data==-ch_num);
    rise_idx = find(eventData.Data== ch_num);
    % some data have rise glitch on first event, if it happen within first 8 evenet, remove it
    if length(fall_idx)+1 == length(rise_idx) && any(rise_idx(1) == (1:8))
        rise_idx = rise_idx(2:end);
    end
    pulse_width = evt_ts(rise_idx) - evt_ts(fall_idx);
    gli_col = zeros(length(pulse_width),1);
    gli_col(pulse_width==1) = 1;
    evt_debug = [evt_debug; [evt_data(fall_idx) evt_ts(fall_idx) evt_ts(rise_idx) pulse_width gli_col]];
end
[~,sort_idx]= sort(evt_debug(:,2)); % sort by timestamp
evt_debug = evt_debug(sort_idx,:); 

% remove glitch
fprintf('%d electrical glitches removed \n',sum(evt_debug(:,4) == 1))
evt_debug(evt_debug(:,4) == 1,:) = [];

% find split events for both rise and fall
evt_fall_ts = evt_debug(:,2);
split_ts_idx = find(diff(evt_fall_ts) == 1);
for i_ts = 1:length(split_ts_idx)
    all_idx = evt_fall_ts(split_ts_idx(i_ts)) == evt_fall_ts | evt_fall_ts(split_ts_idx(i_ts))+1 == evt_fall_ts;
    evt_debug(all_idx,2) = evt_fall_ts(split_ts_idx(i_ts));
end
evt_rise_ts = evt_debug(:,3);
split_ts_idx = find(diff(evt_rise_ts) == 1);
for i_ts = 1:length(split_ts_idx)
    evt_split_ts = evt_rise_ts(split_ts_idx(i_ts));
    all_idx = find(evt_split_ts == evt_rise_ts | evt_split_ts+1 == evt_rise_ts);
    if numel(unique(evt_debug(all_idx,2))) == 1
        evt_debug(all_idx,3) = evt_split_ts;
    else
        disp('e')
    end
end

all_ch_evt = evt_debug;
all_ch_evt(:,6) = all_ch_evt(:,1);
% 1 3 4 5 6 - > 1 2 3 4 5
all_ch_evt(all_ch_evt(:,1) == -1,1) = -2; % replace channel 1 to 2
all_ch_evt(:,1) = abs(all_ch_evt(:,1))-1;  % channel range from 2-6 to 1-5  
uni_ts = unique(all_ch_evt(:,2:3));

% track event channel status
eventType = '00000';
sample_rate = 30000;
event_tt = zeros(length(uni_ts),4);
for i_evt = 1:length(uni_ts)
    ts_on_idx = all_ch_evt(uni_ts(i_evt) == all_ch_evt(:,2),1);
    ts_off_idx = all_ch_evt(uni_ts(i_evt) == all_ch_evt(:,3),1);
    eventType(ts_on_idx) = '1';
    eventType(ts_off_idx) = '0';
    event_tt(i_evt,1) = double(uni_ts(i_evt))/sample_rate;
    event_tt(i_evt,2) = bin2dec(eventType);
    event_tt(i_evt,3) = uni_ts(i_evt);
    pulse_width = all_ch_evt(uni_ts(i_evt) == all_ch_evt(:,2),4);
    if ~isempty(pulse_width)
        event_tt(i_evt,4) = pulse_width(1);
    end
end    

events = struct;
uni_type = unique(event_tt(:,2));
uni_type = uni_type(uni_type>0);
for i_type = 1:length(uni_type)
    events.(sprintf('evt%d',uni_type(i_type))).ts = event_tt(uni_type(i_type) == event_tt(:,2),1);
    if any(event_type_dec == uni_type(i_type))
        events.(sprintf('evt%d',uni_type(i_type))).type = event_record_type{event_type_dec == uni_type(i_type),2};
    else
        events.(sprintf('evt%d',uni_type(i_type))).type ='No defined event type';
        warning(sprintf('no defined type - evt%d',uni_type(i_type)))
    end
end

led_ttl = evt_ts(eventData.Data == 2);
if ~isempty(led_ttl)
    events.evt32.ts = double(led_ttl)./sample_rate;
    events.evt32.type = 'LED TTL';
end

eventnames = fieldnames(events);
for i_type = 1:length(eventnames)
    fprintf('%s: %s - %d events\n', eventnames{i_type}, events.(eventnames{i_type}).type,length(events.(eventnames{i_type}).ts))
end

if nargin == 1
    return
end


% check against mpc data
ani_id = regexp(name,'^\w+(?=_\d{4})','match','once');
date_str = regexp(name,'(?<=_)\d{4}-\d{2}-\d{2}(?=_\d{2})','match','once');

date_str_conv = string(datetime(date_str) , 'MM/dd/yy');
match_mpc_idx = find(strcmp(date_str_conv,{mpcParsed.StartDate}) & strcmp(ani_id,{mpcParsed.Subject}));


if isempty(match_mpc_idx)
    fprintf('NO MATCHING MPC FOUND! \n')
    return
end
mpc_parsed_match = mpcParsed(match_mpc_idx(1));
max_dev_set = 0.3; % maximum allowed deviation
event_check = { 2, 'H'; 17, 'S'; 9, 'W'; 7, 'P'; 11, 'Q'; 19, 'N'; 15, 'O'; 6, 'G';};
max_dev_c = [];
for i_evt = 1:length(event_check)
    evt_c = event_check(i_evt,:);
    evt_oe_ts = events.(sprintf('evt%d',evt_c{1})).ts;
    if evt_c{1} == 2 % trial start
        evt_oe_ts = sort([events.evt22.ts; events.evt23.ts]);
    end    
    evt_mpc_ts = mpc_parsed_match(1).(evt_c{2});
    if length(evt_oe_ts) == length(evt_mpc_ts)
        max_dev = max(abs(diff(evt_oe_ts)-diff(evt_mpc_ts)));
        max_dev_c(i_evt) = max_dev;
        if max_dev < max_dev_set
        else
            fprintf('oe evt %d ts NOT match with mpc %s evt at max deviation %0.2fms %s \n',evt_c{1},evt_c{2},max_dev*1000, events.(sprintf('evt%d',evt_c{1})).type)
        end
    else
        % oe evt miss events
        if numel(evt_oe_ts) < numel(evt_mpc_ts)
            min_num = min([numel(evt_oe_ts) numel(evt_mpc_ts)]);
            max_dev = abs(diff(evt_oe_ts(1:min_num))-diff(evt_mpc_ts(1:min_num)));
            off_idx = find(max_dev>max_dev_set,1);
            off_ts = evt_mpc_ts(off_idx+1) + (evt_oe_ts(1)- evt_mpc_ts(1));
            evt_oe_ts = sort([evt_oe_ts;off_ts]);
            if length(evt_oe_ts) == length(evt_mpc_ts)
                max_dev = max(abs(diff(evt_oe_ts)-diff(evt_mpc_ts)));
                max_dev_c(i_evt) = max_dev;
                fprintf('oe evt %d ts fixed with mpc %s evt %s \n',evt_c{1},evt_c{2}, events.(sprintf('evt%d',evt_c{1})).type)
            end
         % oe evt extra events    
        elseif numel(evt_oe_ts) > numel(evt_mpc_ts)
            min_num = min([numel(evt_oe_ts) numel(evt_mpc_ts)]);
            max_dev = abs(diff(evt_oe_ts(1:min_num))-diff(evt_mpc_ts(1:min_num)));
            off_idx = find(max_dev>max_dev_set,1);
            if numel(off_idx) > 0             
                evt_oe_ts(off_idx+1) = [];
                events.(sprintf('evt%d',evt_c{1})).ts(off_idx+1) = [];
            else % remove last event ts 
                evt_oe_ts(end) = [];
                events.(sprintf('evt%d',evt_c{1})).ts(end) = [];
            end
            if length(evt_oe_ts) == length(evt_mpc_ts)
                max_dev = max(abs(diff(evt_oe_ts)-diff(evt_mpc_ts)));
                max_dev_c(i_evt) = max_dev;
                fprintf('oe evt %d ts fixed with mpc %s evt %s \n',evt_c{1},evt_c{2}, events.(sprintf('evt%d',evt_c{1})).type)
            end
        end

    end
end
if all(max_dev_c < max_dev_set)
    fprintf('all events match with max deviation %0.2fms \n', max(max_dev_c)*1000)
end




