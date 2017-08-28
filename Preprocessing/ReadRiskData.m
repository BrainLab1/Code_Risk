
% initiation: 01.07.2017; created by Bahareh and Saeed; this is the modified
% version of 'fieldtrip_preproc_par.m' by (Amin et al)
% --------------------------------------------------------------------------------------
% This script reads out .bhv file for each session, creats the structure
% array 'event' which includes trial-by-trial behavioral data and task parameters
% and timings (timings are not aligned to events).
% Then reads ot the LFP data from the .ns2 file, sort it into trials,
% alignes the data to events, runs a simple preprocessing routine (notch
% filter, detrendign, low pass filtering) and stores the behavioral and LFP
% data together on the hard disk


clear
close all;
clc;
dbstop if error

main_folder =  'Z:\'; % changed by Bahareh 12-06-2017 (Z: is the shared folder between all computers)

data_folder         = 'data\'; % 
code_folder         = 'code\';
toolbox_folder      = 'toolbox\';

data_dir            = 'Risk\Original Data_Extracted\';
save_dir            = 'Risk\Primary Filtered LFP Data\';

cue_interval        = [-0.2, 1];   % zero is cue onset
reward_interval     = [-0.6, 0.1]; % zero is reward onset

tmp_folder_list     = dir([main_folder data_folder data_dir]);                                 
iniCfg = [];

%% addd required path
addpath( genpath(['Z:\data\Risk\toolbox\fieldtrip-20160904']) );
addpath( genpath(['Z:\data\Risk\toolbox\NPMK']) );
addpath( genpath(['Z:\data\Risk\toolbox\MonkeyLogic']) );
% addpath( genpath(['C:\Users\BrainLab\Dropbox\Jackie Project\shared code\toolbox\fieldtrip-20160904']) );
% addpath( genpath(['C:\Users\BrainLab\Dropbox\Jackie Project\shared code\toolbox\NPMK']) );
% addpath( genpath(['C:\Users\BrainLab\Dropbox\Jackie Project\shared code\toolbox\MonkeyLogic']) );
% addpath( genpath(['C:\Users\BrainLab\Dropbox\Jackie Project\shared code\Git_Repository']) );

%% General config for all sessions

tmp_cfg                 = [];
tmp_cfg.headerformat    = 'blackrock_nsx';
tmp_cfg.data_dir        = [main_folder data_folder data_dir];
tmp_cfg.trialfun        = 'fieldtrip_trialfun_RiskNeural'; % string with function name, see below (default = 'ft_trialfun_general')
tmp_cfg.trialdef.pre    = 1; % time in seconds; this only indicates absolut value, do not use '-'
tmp_cfg.trialdef.post   = 1; % time in seconds

% creat a configuration for all the sessions; iniCfg is a structure array, each row indicates on session and one alignment;
% if you wanna aligne same session data to two different events you need to allocate two rows of the iniCfg to this session
for i = 1:numel(tmp_folder_list) % for each session
    % read out list of all the files related to this session
    tmp_list = tmp_folder_list(i);  
    
    % here we assume that data were collected only from MacDuff and Mojo
    if numel(tmp_list.name) > 3 && (strcmp(tmp_list.name(1:3),'Moj') || strcmp(tmp_list.name(1:3),'Mac'))
        
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
        
        %% Config for aligning on cue
        tmp_cfg.trialdef.eventtype  = 'cue';
        tmp_cfg.trialdef.interval   = cue_interval;
        iniCfg = [iniCfg; tmp_cfg];
        
        %% Config for aligning on reward
        tmp_cfg.trialdef.eventtype  = 'reward';
        tmp_cfg.trialdef.interval   = reward_interval;
        iniCfg = [iniCfg; tmp_cfg];

    end
end

%% Distributed preprocessing

% num_of_workers = 8;
% distcomp.feature('LocalUseMpiexec', true);
% poolobj = parpool('ll', num_of_workers); %% ll refers to Lashgari's Lab!

% % % % poolobj = parpool('local', 4);

% parfor i = 1: length(iniCfg)
for i = 1:length(iniCfg)
    % if data is already preprocessed, continue
    if exist([main_folder data_folder 'Preprocessed - new\' 'Data_' iniCfg(i).session_dir ' (' iniCfg(i).trialdef.eventtype ').mat'], 'file')
        continue;
    end
    
    
    try
        cfg = ft_definetrial(iniCfg(i));
        data = fieldtrip_preproc_session([main_folder data_folder], cfg);
        save ([main_folder data_folder save_dir 'Data_' iniCfg(i).session_dir ' (' iniCfg(i).trialdef.eventtype ').mat'], 'data')
        clear data cfg
    catch err
        warning([err.message ' in session: ' iniCfg(i).session_dir]);
        continue;
    end
end

% % % % delete(poolobj);