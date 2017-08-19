
% last update: Saeed; some bugs are fixed and outliers are named according
% to the parameters. two outlier indices are added 'OutCueOnOffZ2' and
% 'OutAqToRwdZ2' where Z2 means z-score for the outlier removal is 2
% Saeed; this code removes the outliers of different trial timings
% according to the specified properties by the user. it is remarkable that,
% this code add a new field to the saved Behavioral data to specify which
% trials are outliers '0' or not '1'.

clear
close all
clc
dbstop if error

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
main_folder          = 'Z:\'; 
data_folder          = 'data\'; % 
original_data_folder = 'Risk\Original Data_Extracted\';
save_dir             = 'Z:\data\Risk\Behavior\';

%% Read out list of all the files related to this session
dataPath = [main_folder data_folder original_data_folder];
allFiles = dir(dataPath);    

tmp = [];
for i = 1:length(allFiles)
    if  numel(allFiles(i).name) < 3 
        tmp = [tmp; i];
    end
end
allFiles(tmp)=[];
clear i tmp
%%

if 1
    
Events = [];
for ses = 1:numel(allFiles) % for each session
    % select the right file
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
        
        dir_name{ses} = new_cfg.session_dir;

        % convert the event tructure to table
        eventTable = struct2table(new_cfg.event);
        eventTable.SessionID = (ses)*ones(size(eventTable,1),1);
        Events = [Events;eventTable];
        
        clear eventTable new_cfg

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
% outlierfield.trial_cueonset_to_cueoffset = 1;
% outlierfield.trial_cueoffset_to_targetacquired = 1;
% outlierfield.trial_targeton_to_targetacquired = 1;
% outlierfield.trial_target_to_reward = 1;
% outlierfield.trial_reward_duration = 1;
% outlierfield.trial_reaction_time = 1;

t1 = [Events.TrialErrorCode == 0]; % binary vector of trial indices with successful trials '1' and error trials '0'
indx = find(Events.TrialErrorCode == 0);% the indices of successful trials
c = Events(indx,:);% table of successful trials

% detection of outliers from CueOn to CueOff timing
outlierfield.trial_cueonset_to_cueoffset = 1;
[ ~,outlier_indices ] = Outlier_remove( table2struct(c),outlierfield );
for i = 1:numel(outlier_indices)
    outlier_index(i) = indx(outlier_indices(i));
end
Events.OutCueOnOffZ2 = zeros(size(Events,1),1);
Events.OutCueOnOffZ2(outlier_index) = 1;
clear outlier_indices outlier_index

% detection of outliers from Target Acuired to Reward On timing
outlierfield.trial_target_to_reward = 1;
[ ~,outlier_indices ] = Outlier_remove( table2struct(c),outlierfield );
for i = 1:numel(outlier_indices)
    outlier_index(i) = indx(outlier_indices(i));
end
Events.OutAqToRwdZ2 = zeros(size(Events,1),1);
Events.OutAqToRwdZ2(outlier_index) = 1;
clear outlier_indices outlier_index

for j = 1:length(allFiles)
    trl = find(Events.SessionID == j);
    load ([save_dir 'Bhv_' dir_name{j} '.mat'])
    new_cfg.event = table2struct(Events(trl,:));
    new_cfg.event = rmfield(new_cfg.event,'SessionID');
    
    save ([save_dir 'Bhv_' dir_name{j} '.mat'], 'new_cfg')
    clear trl new_cfg
end

