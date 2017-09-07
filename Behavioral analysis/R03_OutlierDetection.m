
% last update: by Saeed; detection of RT outliers with zscore of log is
% added
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
original_data_folder = 'Risk\Behavior\';
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
    
Events = [];
for ses = 1:numel(allFiles) % for each session
    % select the right file
    sessionFolder = allFiles(ses).name;
    load([main_folder data_folder original_data_folder sessionFolder])
    
    % convert the event tructure to table
    eventTable = struct2table(new_cfg.event);
    eventTable.SessionID = (ses)*ones(size(eventTable,1),1);
    Events = [Events;eventTable];
    
    clear eventTable new_cfg
end
clear ses

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
clear outlier_indices outlier_index outlierfield

% detection of outliers from Target Acuired to Reward On timing
outlierfield.trial_target_to_reward = 1;
[ ~,outlier_indices ] = Outlier_remove( table2struct(c),outlierfield );
for i = 1:numel(outlier_indices)
    outlier_index(i) = indx(outlier_indices(i));
end
Events.OutAqToRwdZ2 = zeros(size(Events,1),1);
Events.OutAqToRwdZ2(outlier_index) = 1;
clear outlier_indices outlier_index outlierfield

% detection of outliers from Reaction Time
outlierfield.trial_reaction_time = 1;
[ ~,outlier_indices ] = Outlier_remove( table2struct(c),outlierfield );
for i = 1:numel(outlier_indices)
    outlier_index(i) = indx(outlier_indices(i));
end
Events.OutRTZ2 = zeros(size(Events,1),1);
Events.OutRTZ2(outlier_index) = 1;
clear outlier_indices outlier_index outlierfield


for j = 1:length(allFiles)
    trl = find(Events.SessionID == j);
    load([main_folder data_folder original_data_folder allFiles(j).name])
    new_cfg.event = table2struct(Events(trl,:));
    new_cfg.event = rmfield(new_cfg.event,'SessionID');
    
    save([main_folder data_folder original_data_folder allFiles(j).name],'new_cfg')
    clear trl new_cfg
end

