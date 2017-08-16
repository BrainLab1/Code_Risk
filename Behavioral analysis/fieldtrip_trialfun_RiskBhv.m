
% last update 16.08.2017 by Bahareh:  event times and reaction times corrected for display latency, added to the output.
% last update 08.08.2017:by Saeed; all the information about alignment event are
% removed and one new field named 'Photodiode' is added. three fields named
% 'offset', 'sample' and 'duration' are removed as they are extra for behavioral data
% analysis
% last update 26.07.2017: new fields added to the output structure: 'EyeSignal', 'LickSignal', 'PupilSize', 'ActualEventTime'  
% last update 25.07.2017: new fields added to the output structure: 'value' changed to 'ActualRewardValue'; 'RewardVariance' added 
% created by Bahareh, 18.07.2017

function [trl,event] = fieldtrip_trialfun_RiskBhv(cfg)
% Finding trial and event information using behavioral data and header file
% of the blackrock (ns2) file.

% Created by Amin (nejatbakhsh.amin@gmail.com), modifyed by Bahareh for the Risk experiment
    

    hdr = ft_read_header(cfg.headerfile, 'headerformat', 'blackrock_nsx');
    bhv_dir = cfg.data_dir;
    
    trl = [];
    event = [];

    pretrig  = -round(cfg.trialdef.pre  * hdr.Fs);
    posttrig =  round(cfg.trialdef.post * hdr.Fs);
    
    if isfield(cfg.trialdef,'interval')
        interval_begin = round(cfg.trialdef.interval(1) * hdr.Fs);
        interval_end = round(cfg.trialdef.interval(2) * hdr.Fs);
    end
    
    tmp_list = cfg.session_dir;

    bhv_name = dir([bhv_dir tmp_list '\*.bhv']);
    onoff_name = dir([bhv_dir tmp_list '\' '*onoff']);

    if length(bhv_name) > 1
        warning(['Multiple BHV files in folder ' [bhv_dir tmp_list]]);
        return;
    end

    if length(onoff_name) > 1
        warning(['Multiple ONOFF files in folder ' [bhv_dir tmp_list]]);
        return;
    end

    bhv_file        = bhv_read([bhv_dir tmp_list '\' bhv_name.name]);
    onoff_file      = fileread([bhv_dir tmp_list '\' onoff_name.name]);

    tmp_onoff       = regexp(onoff_file, '{(?<on>.*)}.*{(?<off>.*)}', 'names');
    onoff           = {round(cellfun(@str2num, regexp(tmp_onoff.on, '[0-9]+', 'match'))/30), ...
                       round(cellfun(@str2num, regexp(tmp_onoff.off, '[0-9]+', 'match'))/30)};

    if length(onoff{1}) < length(bhv_file.ConditionNumber)
        warning(['The number of elements in ONOFF file is lower than the number of trials in BHV file in folder ' [bhv_dir tmp_list]]);
        return;
    end
    
     %% Trials
    timing_cue      = [];
    timing_reward   = [];
    

    for j = 1: length(bhv_file.ConditionNumber)
        cue         = bhv_file.CodeTimes{j}(bhv_file.CodeNumbers{j} == 19); if isempty(cue), cue = -Inf; end
        start       = bhv_file.CodeTimes{j}(bhv_file.CodeNumbers{j} == 1); if isempty(start), start = -Inf; end
        target_aq   = bhv_file.CodeTimes{j}(bhv_file.CodeNumbers{j} == 14); if isempty(target_aq), target_aq = -Inf; end
        
        timing_cue  = [timing_cue; cue - start];
        timing_reward = [timing_reward; target_aq - start + 377] ;
    end

    if ~isfield(cfg.trialdef,'eventtype') || isempty(cfg.trialdef.eventtype)
        new_trials  = [onoff{1}(1: length(bhv_file.ConditionNumber))' + pretrig, ...
                      onoff{1}(1: length(bhv_file.ConditionNumber))' + posttrig, ...
                      ones(length(bhv_file.ConditionNumber), 1) * (pretrig)];
    elseif strcmp(cfg.trialdef.eventtype, 'cue')
        new_trials  = [onoff{1}(1: length(bhv_file.ConditionNumber))' + timing_cue + interval_begin + pretrig, ...
                      onoff{1}(1: length(bhv_file.ConditionNumber))' + timing_cue + interval_end + posttrig, ...
                      ones(length(bhv_file.ConditionNumber), 1) * (interval_begin + pretrig)];
    elseif strcmp(cfg.trialdef.eventtype, 'reward')
        new_trials  = [onoff{1}(1: length(bhv_file.ConditionNumber))' + timing_reward + interval_begin + pretrig, ...
                      onoff{1}(1: length(bhv_file.ConditionNumber))' + timing_reward + interval_end + posttrig, ...
                      ones(length(bhv_file.ConditionNumber), 1) * (interval_begin + pretrig)];
    end

    trl = [trl; new_trials];
    
    %% Pass actual event times/samples (not aligned) to the output
    actualEventTime = {};
    eventCodes = [9 1 2 6 19 20 21 14 22 100 18]';
    eventNames = [{'MLTrialStart'}; {'TrialGateOn'}; {'FixationOn'}; {'FixationAcquired'}; {'CueOnset'}; {'CueOffset'};...
                  {'TargetOnset'}; {'TargetAcquired'}; {'TargetOffset'}; {'TrialGateOff'}; {'MLTrialEnd'}];
   tr = 1;
   while tr <= length(bhv_file.ConditionNumber)  % for each trial
       ev = 1;
       while ev <= length(eventCodes)
            if ~isempty( bhv_file.CodeTimes{tr}(bhv_file.CodeNumbers{tr} == eventCodes(ev)) )
                eval(['actualEventTime{tr,1}.' eventNames{ev} ' = bhv_file.CodeTimes{tr}(bhv_file.CodeNumbers{tr} == eventCodes(ev));']);
            else 
                eval(['actualEventTime{tr,1}.' eventNames{ev} ' = NaN;']);
            end
            ev = ev+1;
       end  
       tr = tr+1;
   end
   clear tmp tr ev eventNames eventCodes
   
    %% Extract the PhotoDiod, EyeSignal, LickingSignal and PupilSize 
    eyeSignal   = {};
    lickSignal  = {};
    pupilSize   = {};
    photoDiode  = {};
    for tr = 1: length(bhv_file.ConditionNumber)  % for each trial
        eyeSignal{tr,1}  = bhv_file.AnalogData{tr}.EyeSignal;    
        lickSignal{tr,1} = bhv_file.AnalogData{tr}.General.Gen1;
        pupilSize{tr,1}  = struct('Height', bhv_file.AnalogData{tr}.General.Gen2, ...
                                  'Width', bhv_file.AnalogData{tr}.General.Gen3);
        photoDiode{tr,1} = bhv_file.AnalogData{tr}.PhotoDiode;             
    end
    clear tr
    
    %% Extract useful info from the .bhv file and pass it to output
    condition_structs   = [bhv_file.InfoByCond{:}]';
    condition_names     = reshape({condition_structs.Stim2}, [length(condition_structs) 1]);
    conditions          = condition_names(bhv_file.ConditionNumber);  % task condition for each trial

    reward_structs      = [bhv_file.RewardRecord];    
    c                   = squeeze(struct2cell(reward_structs))' ;  % size: cell{numTrials , 2}
    rewards             = cellfun(@(x,y) sum(y-x), c(:,1), c(:,2), 'UniformOutput', 0);
    clear c
    
    % find the expected reward (mean reward) for each task condition using the task conditionname! 
    % use bhv_file.TaskObject(bhv_file.ConditionNumber, 3) which looks like 'pic(F3,8.000,0.000...'
    expected_rewards = num2cell( cellfun(@(x) str2num(x(6)), bhv_file.TaskObject(bhv_file.ConditionNumber, 3)) ); % 3rd column corresponds to Cue porperty

    goods               = ([bhv_file.TrialError] == 0); % successful trials
    pre_goods           =  num2cell([0; goods(1:end-1)]); % success status of the previous trial
    
    % find the indices for all ',' character in the Cue object property; for all trials at once; 
    % Cue object property looks like this:  bhv_file.TaskObject(trialNumber, 3) = pic(F6,8.000,0.000,120.000,30.000)
    [startIndex,~] = regexp( bhv_file.TaskObject(bhv_file.ConditionNumber, 3),','); % size(startIndex) = {NumberOfTrials , 1}; startIndex{trialNumber} = [1 x 4] matrix
    % convert the startIndex into a matrix
    startIndex = cell2mat(startIndex);  % size(startIndex) = [NumberOfTrials x 4]
    % now get the Cue position
    cue_positions       = num2cell( cellfun(@(x,idx) sign(str2num(x( idx(1)+1:idx(2)-1 ))), bhv_file.TaskObject(bhv_file.ConditionNumber, 3), mat2cell(startIndex,ones(1,size(startIndex,1)))) ); % 3rd column corresponds to Cue porperty
    clear startIndex

    % find the indices for all ',' character in the Target object property; for all trials at once; 
    % Target object property looks like this:  bhv_file.TaskObject(trialNumber, 5) = pic(T1,8.000,0.000,50.000,50.000)
    [startIndex,~] = regexp( bhv_file.TaskObject(bhv_file.ConditionNumber, 5),','); % size(startIndex) = {NumberOfTrials , 1}; startIndex{trialNumber} = [1 x 4] matrix
    % convert the startIndex into a matrix
    startIndex = cell2mat(startIndex);  % size(startIndex) = [NumberOfTrials x 4]
    % now get the Target position
    target_positions       = num2cell( cellfun(@(x,idx) sign(str2num(x( idx(1)+1:idx(2)-1 ))), bhv_file.TaskObject(bhv_file.ConditionNumber, 5), mat2cell(startIndex,ones(1,size(startIndex,1)))) ); % 3rd column corresponds to Cue porperty
    clear startIndex

    % find the reward variance for each trial
    rewardVariance = {};
    trIdx = 1;
    while trIdx <= length(conditions)
        switch(conditions{trIdx})
            case {'F3','F6','F9'}
                rewardVariance{trIdx,1} = 0;
            case {'B3L','B6L','B9L'}
                rewardVariance{trIdx,1} = 1;
            case {'B3H','B6H','B9H'}
                rewardVariance{trIdx,1} = 4;
            otherwise
                rewardVariance{trIdx,1} = NaN;
        end
        trIdx = trIdx+1;
    end
    clear trIdx

    reawardTime         = squeeze(struct2cell(reward_structs))';
    reaction_time       = mat2cell(bhv_file.ReactionTime',ones(size(bhv_file.ReactionTime',1),1));
    trial_error_code    = mat2cell(bhv_file.TrialError,ones(size(bhv_file.TrialError,1),1)); 
    
    %% Correct ActualEventTime and ReactionTimes for display latency and pass it to the output
    diodEventTime = {};
    diodeRT = {};  % this will include corrected reaction times
    
    for tr = 1: length(bhv_file.ConditionNumber)  % for each trial
        % get the time bin in which diode signal changes
        diodChangTimBin = find(diff(photoDiode{tr,1})>2.5 | diff(photoDiode{tr,1})<-2.5);

        % treat successful and failed trials differently!
        switch trial_error_code{tr}
            case 0  % for successful trials, diode changes 6 times
                % correct event times
                diodEventTime{tr,1}.FixationOn   = diodChangTimBin(1);
                diodEventTime{tr,1}.CueOnset     = diodChangTimBin(2);
                diodEventTime{tr,1}.CueOffset    = diodChangTimBin(3);
                diodEventTime{tr,1}.TargetOnset  = diodChangTimBin(4);
                % the 5th change of the diode signal occures at the time of reward!?
                diodEventTime{tr,1}.TargetOffset = diodChangTimBin(6);
                % correct reaction time
                diodeRT{tr,1} = reaction_time{tr} - (diodEventTime{tr,1}.TargetOnset - actualEventTime{tr,1}.TargetOnset);
                
            otherwise
                % correcte available event times
                if ~isnan(actualEventTime{tr,1}.FixationOn)
                    diodEventTime{tr,1}.FixationOn   = diodChangTimBin(1);
                    if ~isnan(actualEventTime{tr,1}.CueOnset)
                        diodEventTime{tr,1}.CueOnset     = diodChangTimBin(2);
                        if ~isnan(actualEventTime{tr,1}.CueOffset)
                            diodEventTime{tr,1}.CueOffset    = diodChangTimBin(3);
                            if ~isnan(actualEventTime{tr,1}.TargetOnset)
                                diodEventTime{tr,1}.TargetOnset  = diodChangTimBin(4);
                                if (~isnan(actualEventTime{tr,1}.TargetOffset) & ~isnan(actualEventTime{tr,1}.TargetAcquired))
                                    diodEventTime{tr,1}.TargetOffset = diodChangTimBin(6);
                                end
                            end
                        end
                    end
                end
                % correct reaction time if available
                if isfield(diodEventTime{tr,1}, 'TargetOnset')
                    diodeRT{tr,1} = reaction_time{tr} - (diodEventTime{tr,1}.TargetOnset - actualEventTime{tr,1}.TargetOnset);
                else
                    diodeRT{tr,1} = NaN;
                end
            end
    end
	clear tr
    
    
    %% form the output on 
    new_events          = struct('type', conditions, 'TotalRewardTime', rewards, ...
                            'expected_reward', expected_rewards, 'RewardVariance', rewardVariance, ...
                            'RewardOnTime',  reawardTime(:,1), 'RewardOffTime',  reawardTime(:,2), 'ReactionTime',  reaction_time, 'DiodeReactionTime', diodeRT,...
                            'ActualEventTime', actualEventTime, 'DiodeEventTime', diodEventTime, 'cue_pos', cue_positions, 'target_pos', target_positions, ...
                            'pre_good', pre_goods, 'TrialErrorCode', trial_error_code, ... 
                            'EyeSignal', eyeSignal, 'LickSignal', lickSignal, 'PupilSize', pupilSize, 'PhotoDiode', photoDiode);
    
    event               = [event; new_events];

%     indices             = find([bhv_file.TrialError]' == 0 & ...
%                                [event.duration] > 0 & ...
%                                [event.sample] > 0 & ...
%                                [event.sample] + [event.duration] < hdr.nSamples);
%     event               = event(indices);
%     trl                 = trl(indices, :);
    
end