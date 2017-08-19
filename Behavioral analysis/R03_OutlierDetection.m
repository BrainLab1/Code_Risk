

% Saeed; this code removes the outliers of different trial timings
% according to the specified properties by the user. it is remarkable that,
% this code add a new field to the saved Behavioral data to specify which
% trials are outliers '0' or not '1'.

clear all;
close all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
main_folder          = 'Z:\'; 
data_folder          = 'data\'; % 
original_data_folder = 'Risk\Original Data_Extracted\';
save_dir             = 'Z:\data\Risk\Behavior\';

% Set the options
groupingCriteria       = 'subjectID' ; % this parameter shows based on what criteria trials are grouped.

%% Read out list of all the files related to this session
dataPath = [main_folder data_folder original_data_folder];
allFiles = dir(dataPath);    

%%

if 1
    
Events = [];
allSesSuccessRT = {}; % reaction time of all successful trials from all sessions and the group identifier -> [ {[reatTime saccadeTarget]} , {monkeyName}, {zscore(reactionTime)} ] 
for ses = 1:numel(allFiles) % for each session
    % select the right file
    if numel(allFiles(ses).name) > 3   
        sessionFolder = allFiles(ses).name;
        bhvFileInfo = dir([dataPath sessionFolder '\*.bhv']); 
        bhvFilePath = [dataPath sessionFolder '\' bhvFileInfo.name];
        
        lfpFileInfo    = dir([dataPath sessionFolder '\*.ns2']);   
        if length(lfpFileInfo) > 1
            warning(['Multiple NS files in folder ' [dataPath sessionFolder]]);
            continue;
        end       
        lfpFilePath = [dataPath sessionFolder '\' lfpFileInfo.name];

        tmp_cfg                 = [];
        tmp_cfg.headerformat    = 'blackrock_nsx';
        tmp_cfg.data_dir        = dataPath;
        tmp_cfg.trialfun        = 'fieldtrip_trialfun_RiskBhv'; % string with function name, see below (default = 'ft_trialfun_general')
        tmp_cfg.trialdef.pre    = 1; % time in seconds; this only indicates absolut value, do not use '-'
        tmp_cfg.trialdef.post   = 1; % time in seconds
        
        % Specific config for one session
        tmp_cfg.headerfile  = lfpFilePath;
        tmp_cfg.dataset     = lfpFilePath;
        tmp_cfg.session_dir = sessionFolder;
        
        tmp_cfg.trialdef.eventtype  = 'cue';
        tmp_cfg.trialdef.interval   = [-0.1 1];

        % read ou the data the same way it is done in ReadRiskDataBhv.m
        new_cfg = ft_definetrial(tmp_cfg);
        
        dir_name{ses-2} = new_cfg.session_dir;
        
        % group trials for 18 task conditions
%         [output] = GroupTrials(new_cfg.event, groupingCriteria);

        % convert the event tructure to table
        eventTable = struct2table(new_cfg.event);
        eventTable.SessionID = (ses-2)*ones(size(eventTable,1),1);
        Events = [Events;eventTable];
        
        clear eventTable new_cfg

    end
end
clear ses


end


%% specify the parameter which you would like to remove outliers from it

% outlierfield.trial_start_time = 1;
% outlierfield.trial_fixation_on = 1;
% outlierfield.trial_fixation_Acquired = 1;
% outlierfield.trial_cue_onset = 1;
% outlierfield.trial_cue_offset = 1;
% outlierfield.trial_target_onset = 1;
% outlierfield.trial_target_acquired = 1;
% outlierfield.trial_target_offset = 1;
% outlierfield.trial_end_time = 1;
% outlierfield.trial_reward_onset = 1;
% outlierfield.trial_reward_offset = 1;
% outlierfield.trial_trial_length = 1;
% outlierfield.trial_fixation_to_acquired = 1;
outlierfield.trial_cueonset_to_cueoffset = 1;
% outlierfield.trial_cueoffset_to_targetacquired = 1;
% outlierfield.trial_targeton_to_targetacquired = 1;
outlierfield.trial_target_to_reward = 1;
% outlierfield.trial_reward_duration = 1;
% outlierfield.trial_reaction_time = 1;

t1 = [Events.TrialErrorCode == 0];
t2 = find(t1 > 0);

indx = find(Events.TrialErrorCode == 0);
c = Events(indx,:);

[ totalTrials_new,outlier_indices ] = Outlier_remove( table2struct(c),outlierfield );


for i = 1:numel(outlier_indices)
    outlier_index(i) = t2(outlier_indices(i));
    ses_id(i) = Events.SessionID(outlier_index(i));
end

num_session = numel(unique([totalTrials_new.SessionID]));

for j = 1:num_session
    load ([save_dir 'Bhv_' dir_name{j} '.mat'])
    Q = struct2table(new_cfg.event);
    Q.OutlierIndex = zeros(size(new_cfg.event,1),1);
    
    q1 = find(ses_id == j);
    for i = 1:numel(q1)
        Q.OutlierIndex(outlier_index(q1(i))) = 1;
    end
    
    new_cfg.event = table2struct(Q);
    
    save ([save_dir 'Bhv_' dir_name{j} '.mat'], 'new_cfg')
    clear new_cfg Q q1
end

