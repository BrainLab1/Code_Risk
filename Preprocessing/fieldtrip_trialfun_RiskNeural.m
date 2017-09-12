
% created by Bahareh and Saeed, 28.08.2017
% unlike fieldtrip_trialfun_RiskBhv which doesn't align the data/timings
% and sends all the trials to the output (regardless of success state),
% this function only passes successful trials to and also alignes all the
% timings and continues behavioral data


function [trl,event] = fieldtrip_trialfun_RiskNeural(cfg)
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
    cue_time = [];
    reward_time = [];

    for j = 1: length(bhv_file.ConditionNumber)
        cue         = bhv_file.CodeTimes{j}(bhv_file.CodeNumbers{j} == 19); if isempty(cue), cue = -Inf; end
        start       = bhv_file.CodeTimes{j}(bhv_file.CodeNumbers{j} == 1); if isempty(start), start = -Inf; end
        target_aq   = bhv_file.CodeTimes{j}(bhv_file.CodeNumbers{j} == 14); if isempty(target_aq), target_aq = -Inf; end
        
        timing_cue  = [timing_cue; cue - start];
        cue_time = [cue_time; cue];
        timing_reward = [timing_reward; target_aq - start + 377] ;
        reward_time = [reward_time; target_aq + 377];
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
    
    %% Pass actual event times/samples (not aligned) and aligned actual event time to the output; 
    alignedEventTime = {};
    actualEventTime = {};
    eventCodes = [9 1 2 6 19 20 21 14 22 100 18]';
    eventNames = [{'MLTrialStart'}; {'TrialGateOn'}; {'FixationOn'}; {'FixationAcquired'}; {'CueOnset'}; {'CueOffset'};...
                  {'TargetOnset'}; {'TargetAcquired'}; {'TargetOffset'}; {'TrialGateOff'}; {'MLTrialEnd'}];
   tr = 1;
   if strcmp(cfg.trialdef.eventtype, 'cue')
       while tr <= length(bhv_file.ConditionNumber)  % for each trial
           ev = 1;
           while ev <= length(eventCodes)
               if ~isempty( bhv_file.CodeTimes{tr}(bhv_file.CodeNumbers{tr} == eventCodes(ev)) )
                   eval(['alignedEventTime{tr,1}.' eventNames{ev} ...
                       ' = bhv_file.CodeTimes{tr}(bhv_file.CodeNumbers{tr} == eventCodes(ev)) - cue_time(tr);']);
                   eval(['actualEventTime{tr,1}.' eventNames{ev} ...
                       ' = bhv_file.CodeTimes{tr}(bhv_file.CodeNumbers{tr} == eventCodes(ev));']);
               else
                   eval(['alignedEventTime{tr,1}.' eventNames{ev} ' = NaN;']);
                   eval(['actualEventTime{tr,1}.' eventNames{ev} ' = NaN;']);
               end
               ev = ev+1;
           end
           tr = tr+1;
       end
       clear tmp tr ev eventNames eventCodes
   elseif strcmp(cfg.trialdef.eventtype, 'reward')
       while tr <= length(bhv_file.ConditionNumber)  % for each trial
           ev = 1;
           while ev <= length(eventCodes)
               if ~isempty( bhv_file.CodeTimes{tr}(bhv_file.CodeNumbers{tr} == eventCodes(ev)) )
                   eval(['alignedEventTime{tr,1}.' eventNames{ev} ...
                       ' = bhv_file.CodeTimes{tr}(bhv_file.CodeNumbers{tr} == eventCodes(ev)) - reward_time(tr);']);
                   eval(['actualEventTime{tr,1}.' eventNames{ev} ...
                       ' = bhv_file.CodeTimes{tr}(bhv_file.CodeNumbers{tr} == eventCodes(ev));']);                   
               else
                   eval(['alignedEventTime{tr,1}.' eventNames{ev} ' = NaN;']);
                   eval(['actualEventTime{tr,1}.' eventNames{ev} ' = NaN;']);
               end
               ev = ev+1;
           end
           tr = tr+1;
       end
       clear tmp tr ev eventNames eventCodes
   end
   clear cue_time
   
    %% Assign subject ID for the two monkeys
   if strcmp(bhv_file.SubjectName,'Mojo')
       subjectID = num2cell(ones(size(trl,1),1));
   elseif strcmp(bhv_file.SubjectName,'MacDuff')
       subjectID = num2cell(2*ones(size(trl,1),1));
   end
   
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
    
    %% Correct ActualEventTime and ReactionTimes for display latency and pass it to the output; also aligne the corrected timings and pass it to the output
    diodEventTime = {}; % corrected eventtimes for display latency (not aligned)
    diodeRT = {};  % this will include corrected reaction times
    
    for tr = 1: length(bhv_file.ConditionNumber)  % for each trial
        % get the time bin in which diode signal changes
        diodChangTimBin = find(diff(photoDiode{tr,1})>2.5 | diff(photoDiode{tr,1})<-2.5);
        
        % initialize the diodEventTime for current trial
        diodEventTime{tr,1}.CueOnset = NaN;
        
        if strcmp(cfg.trialdef.eventtype, 'cue')
            if ~isnan(actualEventTime{tr}.CueOnset)  % if CueOnset event occured
                % get the min-distance diode change to the CueOnset
                [~,min_indx] = min(abs(actualEventTime{tr}.CueOnset - diodChangTimBin));
                % take the diode change time as the correct event time for CueOnset
                diodEventTime{tr,1}.CueOnset     = diodChangTimBin(min_indx);
                alignedDiodEventTime{tr,1}.CueOnset = 0;
                
                % correct and aligne FixationOn
                if ~isnan(actualEventTime{tr}.FixationOn)
                    [~,min_indx] = min(abs(actualEventTime{tr}.FixationOn - diodChangTimBin));
                    diodEventTime{tr,1}.FixationOn = diodChangTimBin(min_indx);
                    alignedDiodEventTime{tr,1}.FixationOn = diodChangTimBin(min_indx) - diodEventTime{tr,1}.CueOnset;
                end
                % correct and aligne CueOffset
                if ~isnan(actualEventTime{tr}.CueOffset)
                    [~,min_indx] = min(abs(actualEventTime{tr}.CueOffset - diodChangTimBin));
                    diodEventTime{tr,1}.CueOffset    = diodChangTimBin(min_indx);
                    alignedDiodEventTime{tr,1}.CueOffset = diodChangTimBin(min_indx) - diodEventTime{tr,1}.CueOnset;
                end
                % correct and aligne TargetOnset
                if ~isnan(actualEventTime{tr}.TargetOnset)
                    [~,min_indx] = min(abs(actualEventTime{tr}.TargetOnset - diodChangTimBin));
                    diodEventTime{tr,1}.TargetOnset  = diodChangTimBin(min_indx);
                    alignedDiodEventTime{tr,1}.TargetOnset = diodChangTimBin(min_indx) - diodEventTime{tr,1}.CueOnset;
                end
                % correct and aligne TargetOffset
                if ~isnan(actualEventTime{tr}.TargetOffset)
                    [~,min_indx] = min(abs(actualEventTime{tr}.TargetOffset(1) - diodChangTimBin));
                    diodEventTime{tr,1}.TargetOffset = diodChangTimBin(min_indx);
                    alignedDiodEventTime{tr,1}.TargetOffset = diodChangTimBin(min_indx) - diodEventTime{tr,1}.CueOnset;
                end
            end            

            
        elseif strcmp(cfg.trialdef.eventtype, 'reward')     
            if ~trial_error_code{tr} % if the current trial was rewarded
                % get the min-distance diode change to the FixationOn
                [~,min_indx] = min(abs(actualEventTime{tr}.FixationOn - diodChangTimBin));
                % take the diode change time as the correct event time for FixationOn
                diodEventTime{tr,1}.FixationOn = diodChangTimBin(min_indx);
                alignedDiodEventTime{tr,1}.FixationOn = diodChangTimBin(min_indx) - reward_time(tr);
                
                % correct and aligne CueOnse
                [~,min_indx] = min(abs(actualEventTime{tr}.CueOnset - diodChangTimBin));
                diodEventTime{tr,1}.CueOnset     = diodChangTimBin(min_indx);
                alignedDiodEventTime{tr,1}.CueOnset = diodChangTimBin(min_indx) - reward_time(tr);

                % correct and aligne CueOffset
                [~,min_indx] = min(abs(actualEventTime{tr}.CueOffset - diodChangTimBin));
                diodEventTime{tr,1}.CueOffset    = diodChangTimBin(min_indx);
                alignedDiodEventTime{tr,1}.CueOffset = diodChangTimBin(min_indx) - reward_time(tr);

                % correct and aligne TargetOnset
                [~,min_indx] = min(abs(actualEventTime{tr}.TargetOnset - diodChangTimBin));
                diodEventTime{tr,1}.TargetOnset  = diodChangTimBin(min_indx);
                alignedDiodEventTime{tr,1}.TargetOnset = diodChangTimBin(min_indx) - reward_time(tr);
                
                % correct and aligne TargetOffset                
                [~,min_indx] = min(abs(actualEventTime{tr}.TargetOffset(1) - diodChangTimBin));
                diodEventTime{tr,1}.TargetOffset = diodChangTimBin(min_indx);
                alignedDiodEventTime{tr,1}.TargetOffset = diodChangTimBin(min_indx) - reward_time(tr);
            end
            
        end
              
        % correct reaction time if available
        if isfield(diodEventTime{tr,1}, 'TargetOnset')
            diodeRT{tr,1} = reaction_time{tr} - (diodEventTime{tr,1}.TargetOnset - actualEventTime{tr,1}.TargetOnset);
        else
            diodeRT{tr,1} = NaN;
        end
    end
	clear tr

%% actual trial index
actualIndx = num2cell((1:numel(subjectID))');
num2cell(ones(size(trl,1),1));

%% pre-trial information
    preTrial{1,1}.type            = NaN;
    preTrial{1,1}.ActualTrialIdx  = NaN;
    preTrial{1,1}.TrialErrorCode  = NaN;
    preTrial{1,1}.expected_reward = NaN;
    preTrial{1,1}.RewardVariance  = NaN;
    preTrial{1,1}.TotalRewardTime = NaN;
    preTrial{1,1}.cue_pos         = NaN;
    preTrial{1,1}.target_pos      = NaN;

    for tr = 2: length(bhv_file.ConditionNumber)  % for each trial
    preTrial{tr,1}.type            = conditions{tr-1};
    preTrial{tr,1}.ActualTrialIdx  = actualIndx{tr-1};
    preTrial{tr,1}.TrialErrorCode  = trial_error_code{tr-1};
    preTrial{tr,1}.expected_reward = expected_rewards{tr-1};
    preTrial{tr,1}.RewardVariance  = rewardVariance{tr-1};
    preTrial{tr,1}.TotalRewardTime = rewards{tr-1};
    preTrial{tr,1}.cue_pos         = cue_positions{tr-1};
    preTrial{tr,1}.target_pos      = target_positions{tr-1};
end
    
    %% form the output on 
    new_events          = struct('ActualTrialIdx',actualIndx,'subjectID',subjectID,'type', conditions, 'TotalRewardTime', rewards, ...
                            'expected_reward', expected_rewards, 'RewardVariance', rewardVariance, ...
                            'RewardOnTime',  reawardTime(:,1), 'RewardOffTime',  reawardTime(:,2), 'ReactionTime',  reaction_time, 'DiodeReactionTime', diodeRT,...
                            'ActualEventTime',actualEventTime,'AlignedEventTime', alignedEventTime, 'DiodeEventTime', diodEventTime, 'AlignedDiodeEventTime', alignedDiodEventTime,...
                            'cue_pos', cue_positions, 'target_pos', target_positions, ...
                            'TrialErrorCode', trial_error_code, 'pre_good', pre_goods, 'PreTrial',preTrial, ... 
                            'EyeSignal', eyeSignal, 'LickSignal', lickSignal, 'PupilSize', pupilSize, 'PhotoDiode', photoDiode);
    
    event               = [event; new_events];

    indices = find([bhv_file.TrialError]' == 0);
    
%     indices             = find([bhv_file.TrialError]' == 0 & ...
%                                [event.duration] > 0 & ...
%                                [event.sample] > 0 & ...
%                                [event.sample] + [event.duration] < hdr.nSamples);
    event               = event(indices);
    trl                 = trl(indices, :);
    
end