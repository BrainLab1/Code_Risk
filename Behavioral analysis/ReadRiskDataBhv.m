
% initiation: 18.07.2017; created by Bahare; this scripts reads out the
% behavioral data from the original data files and sorts everything into
% trials. you need to access the original .bhv file. this script calls
% 'fieldtrip_trialfun_RiskBhv.m'. pay attention to the commented lines at
% the very end of the 'fieldtrip_trialfun_RiskBhv.m'
% last update: 08.08.2017; by Saeed;as this function is for extraction of
% behavioral data, information about 'Cue' or 'Reward' alignment are
% removed


clear
close all
clc;
dbstop if error

main_folder =  'Z:\'; % changed by Bahareh 12-06-2017 (Z: is the shared folder between all computers)

data_folder         = 'data\'; % 
code_folder         = 'code\';
toolbox_folder      = 'toolbox\';

data_dir            = 'Risk\Original Data_Extracted\';
save_dir            = 'Risk\Behavior\';

tmp_folder_list     = dir([main_folder data_folder data_dir]);
% % % % load tmp_folder_list
cfg = [];

%% General config for all sessions

tmp_cfg                 = [];
tmp_cfg.headerformat    = 'blackrock_nsx';
tmp_cfg.data_dir        = [main_folder data_folder data_dir];
tmp_cfg.trialfun        = 'fieldtrip_trialfun_RiskBhv'; % string with function name, see below (default = 'ft_trialfun_general')
tmp_cfg.trialdef.pre    = 1; % time in seconds; this only indicates absolut value, do not use '-'
tmp_cfg.trialdef.post   = 1; % time in seconds

% creat a cfg for all the sessions; cfg is a structure array, each row indicates on session and one alignment;
% if you wanna aligne same session data to two different events you need to allocate two rows of the cfg to this session
for i = 1:numel(tmp_folder_list) % for each session
    % read out list of all the files related to this session
    tmp_list = tmp_folder_list(i);  
    
    % here we assume that data were collected only from MacDuff and Mojo
    if numel(tmp_list.name) > 3 && (strcmp(tmp_list.name(1:4),'Mojo') || strcmp(tmp_list.name(1:7),'MacDuff'))
        
        session_dir = tmp_list.name;
        lfp_file    = dir([main_folder data_folder data_dir tmp_list.name '\' '*.ns2']);
        
        if length(lfp_file) > 1
            warning(['Multiple NS files in folder ' [main_folder data_folder data_dir tmp_list.name]]);
            continue;
        end
        
        filepath = [main_folder data_folder data_dir session_dir '\' lfp_file.name];
        
        %% Specific config for one session
        tmp_cfg.headerfile  = filepath;
        tmp_cfg.dataset     = filepath;
        tmp_cfg.session_dir = session_dir;

        cfg = [cfg; tmp_cfg];

    end
end

%% Distributed preprocessing

% num_of_workers = 8;
% distcomp.feature('LocalUseMpiexec', true);
% poolobj = parpool('ll', num_of_workers); %% ll refers to Lashgari's Lab!

% % % % poolobj = parpool('local', 4);

% parfor i = 1: length(cfg)
for i = 1:length(cfg)
    % if data is already preprocessed, continue
    if exist([main_folder data_folder 'Preprocessed - new\' 'Data_' cfg(i).session_dir '.mat'], 'file')
        continue;
    end
    
    try
        new_cfg = ft_definetrial(cfg(i));
        save ([main_folder data_folder save_dir 'Bhv_' cfg(i).session_dir '.mat'], 'new_cfg')
        clear data new_cfg
    catch err
        warning([err.message ' in session: ' cfg(i).session_dir]);
        continue;
    end
end

% % % % delete(poolobj);