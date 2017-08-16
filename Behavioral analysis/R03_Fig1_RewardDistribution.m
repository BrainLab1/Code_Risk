
clear all;
% close all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
main_folder          = 'Z:\'; 
data_folder          = 'data\'; % 
original_data_folder = 'Risk\Original Data_Extracted\';
save_dir             = 'Z:\data\Risk\Behavior\';

% Set the options
groupingCriteria       = 'type' ; % this parameter shows based on what criteria trials are grouped.

%% Read out list of all the files related to this session
dataPath = [main_folder data_folder original_data_folder];
allFiles = dir(dataPath);                                 

%%
% load matlab_fig1
if 1
    
    
numSuccess = [];
RewardTime = cell(1,9);
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
        
        % group trials for 18 task conditions
        [output] = GroupTrials(new_cfg.event, groupingCriteria);
        
        % convert the event tructure to table
        eventTable = struct2table(new_cfg.event);
        
        % get the trial index for successful trials per condition
        for gr = 1:length(output)
            Names{gr} = output(gr).Value;
            output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
            
            tmp = eventTable([output(gr).SuccessTrlIdx],:);
            RewardTime{gr} = [RewardTime{gr},tmp.TotalRewardTime'];
        end
        clear gr output eventTable new_cfg tmp_cfg tmp
    end
end
clear ses


end

%% set where each condition is ploted
names_ind = [7 8 9
             2 4 6
             1 3 5];
       
%% plot each condition in its predefined position
numRows = 3;
numColumns = 3;
figure('Name','Reward Duration')

for rw = 1:numRows
    for cl = 1:numColumns
        % find the frame number in the subplots
        plotIdx = cl+((rw-1)*numColumns);
        subplot( numRows, numColumns,plotIdx ), hold on, box on
        histogram(RewardTime{names_ind(rw,cl)})
        xlim([0,450]),ylim([0,850])
        xlabel('Total Reward Duration')
        switch plotIdx
            case 1
                title('Expected Reward = 3'),ylabel({'Variance = 0';'';'number of trials'})
            case 2
                title('Expected Reward = 6')
            case 3
                title('Expected Reward = 9')
            case 4
                ylabel({'Variance = 1';'';'number of trials'})
            case 7
                ylabel({'Variance = 4';'';'number of trials'})
        end
    end
end


