
% last update: 04.03.2018: previously we were told that data.hdr.TimeStampPerSample was related to sampling rate of the spike waveform
%       and was set to 30*5 as was for the uncertainty data set. it caused wrong trila-sorting of spikes. data.hdr.TimeStampPerSample is 
%       in fact deermining the difference between LFP and Spike sampling rate. in the curren verion of the Risk data LFP sampling rate 
%       is 1K and spike time sampling rate is 30K, so we should set the data.hdr.TimeStampPerSample to 30. The result of spike trial sorting
%       was the same as when we tried it manually!
% last update: 07.02.2018
% This script 1- loads the LFP data for each session, 2- reads out the sorted spike data for each electrode on the session, 
% cuts the spike train into trials and aligns spike time to events using data.cfg.trl. if the event was Cue,
% first the CueOn time will be corrected using the diode rise time. 

% -- ATTENTION (Bahareh!) --
% pay attention to the lines with following commands. on this line
% asignments are different for different tasks and must be set manually! 
% 
% data.hdr.Fs = ...

%%

clear all;
close all;
clc;
dbstop if error 

main_folder = 'D:\from Z\'; % changed by Bahareh 12-06-2017 (Z: is the shared folder between all computers)

data_folder         = 'data\';
code_folder         = 'code\';
toolbox_folder      = 'toolbox\';

lfp_dir             = 'Z:\data\Risk\Clean Data\';
spike_dir           = 'X:\data\Original Data_Extracted\';
save_dir            = 'X:\data\Cut and Aligned Spk\';

data_type = {'(cue)'};  % events which were used for previouse alignment 

% add toolboxes here
addpath([main_folder code_folder]);

% add the required toolboxes to the MATLAB active path
% addpath([main_folder toolbox_folder 'fieldtrip-20160904']); ft_defaults;
% addpath([main_folder toolbox_folder 'NPMK']); installNPMK;
% addpath([main_folder toolbox_folder 'MonkeyLogic']);

%% loop over data types (cue, reward)
for DT = 1:numel(data_type) % for each alignment
tmp_folder_list = dir([lfp_dir '*' data_type{DT} '.mat']); % Bahareh: this is the list of all field data files

%% loop over sessions
for i = 12:numel(tmp_folder_list) % for each session
    %% load spike and cleaned (200Hz for uncertainty task; 1K for risk task) lfp data like field trip
    load([lfp_dir tmp_folder_list(i).name])  % Bahareh: this is the selected LFP .mat file

    % read out the monkey name
    switch tmp_folder_list(i).name(6:8)
        case 'Mac'
            Monkey = 'MacDuff';
        case 'Moj'
            Monkey = 'Mojo';
    end
        
    spike = []; % Bahareh: define the new 'spike' structure
    spike.hdr = [];
    spike.hdr.FileHeader.Frequency = 30000;
    
    spike_counter = 0; % Bahareh: this is counting the selected sorted units on one electrode 
    spike_list = dir([spike_dir tmp_folder_list(i).name(6:end-5-numel(data_type{DT})) '\*.mat']); % Bahareh: this is list of all the spike .mat files corresponding to the selected field file!
    
    % if there is no isolated unit on this session, continue the loop and go to the next session
    if isempty(spike_list)
        continue;
    end
    
    for SL = 1:numel(spike_list) % for each electrode (each spike file contains all the spikes on one electrode)
        % Bahare: there might be other .mat files under this folder. only load the spike related .mat files.  
        if ~strcmp(spike_list(SL).name(1:10), 'SortSpikes')
            continue;
        end
        
        % load variable Expression1 for this electrode
        load([spike_dir tmp_folder_list(i).name(6:end-5-numel(data_type{DT})) '\' spike_list(SL).name]); % Bahareh: this is one of the spike files in the list that corresponds to the selected LFP file
        % only take the selected spikes (selected during sorting)
        Expression1 = Expression1(Expression1(:,1)==1, :); % Bahareh: each spike file contains a matrix named Expression1, each row represents one spike
        
        group_spikes = Expression1;
        channel_ID = ElecNo2ChanNo_Risk(Monkey, unique(Expression1(:,4)));
        
        groups = unique(Expression1(:,2)); % Bahareh: this indicates number of units on each electrode/channel 
        for GN = 1:numel(groups)  % for each sorted unit
            spike_counter = spike_counter + 1; 
            group_spikes = Expression1(Expression1(:,2)==groups(GN), :);
            spike.label{spike_counter} = [data.label{channel_ID} '-spk-' num2str(GN) '-ch-' num2str(channel_ID)]; % ['elec'array-elec on the array-'spk'-unit number-'ch'-channel/pin number] ??????????
            spike.timestamp{spike_counter} = group_spikes(:,3)';
            spike.waveform{spike_counter} = reshape(group_spikes(:,5:end)',1,size(group_spikes(:,5:end),2),size(group_spikes,1));
            spike.unit{spike_counter} = spike.timestamp{spike_counter}*NaN;  % this is just initialization. true values will be added via 'ft_spike_select(cfg, spike)', see bellow.
        end
        clear Expression1
    end
    
    %% append spike to LFP data
    % change the spike data into the fieldtrip format
    cfg              = [];
    cfg.spikechannel = spike.label; % this is the unit name
    spike            = ft_spike_select(cfg, spike);
    
    % add some more field to the LFP .hdr field, and modify some fields (here due to resampling)
    data.hdr.FirstTimeStamp = 0;
    data.hdr.TimeStampPerSample = 30; % waveform sampling freq.
    
    data.cfg.trl = round(recursive_read(data, 'trl')/(data.hdr.Fs/data.fsample) ); % this will be used bellow for cutting spike data into trials
    data.hdr.Fs = 1000;  % 1000 for risk data; 200 for uncertainty data
    data.hdr.nSamples = round(data.hdr.nSamples);
    
    % turn the spike data into trials
    cfg           = []; 
    cfg.hdr       = data.hdr; % contains information for conversion of samples to timestamps
    cfg.trlunit   = 'samples';
    cfg.trl       = data.cfg.trl; % now in samples
    
    % ---- apply the display latency to the Cue event time for aligning the spike times
    if strcmp(data_type(DT),'(cue)')
        % read out the event structure from the lfp data file
        event = recursive_read(data,'event');
        % calculate diplay latency for Cue Onset event in all trials
        displayLatency = cellfun(@(x,y)  x.CueOnset - y.CueOnset  , {event.DiodeEventTime}' , {event.ActualEventTime}');
        % correct the start/end time of the trial (to cut the spike trains) for the display latency
        cfg.trl(:,1:2) = cfg.trl(:,1:2) + displayLatency;
    end
    
    % cut the spike train into trials and align spike times to the Cue  
    spikeTrials   = ft_spike_maketrials(cfg,spike); % number of spikes in spikeTrials is less tha Expression1 because between trial spikes were removed
    % $$ ATTENTION $$ spikeTrials.trial{un} should indicate each spike of the unit happend in which trial, how ever this may not be the actual trial number because
    % before getting this parameter, we had selected the correct trials. this just the trial index within the correct trials 
    
    spikeTrials.cfg.event = recursive_read(data, 'event');
    
%     data_lfp = data;
%     data_all = ft_appendspike([],data, spikeTrials);

    
    %% save data here
    save([save_dir tmp_folder_list(i).name], 'spikeTrials')
    clear data spikeTrials spike
    
end  % end of sessions...
end  % end of alignments

