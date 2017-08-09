

clear all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
main_folder          = 'Z:\'; 
data_folder          = 'data\'; % 
original_data_folder = 'Risk\Original Data_Extracted\';
save_dir             = 'Z:\data\Risk\Behavior\';

colorCode = [  204 204 255  % F3
               255 153 153  % F6
               153 255 204   % F9
               102 102 255   % B3L
               255 51  51    % B6L
               51  255 102   % B9L              
               0   0   255   % B3H
               204 0   0     % B6H
               0  153 153]/255; % B9H   
           
rtOutlierFactor = 2; % outliers are trials out of this intervel [-1 1]*trOutlierFactor*Std; Std for each monkey is different
lickBestTime = [-0.35 -0.2]; % this is the time (relative to event) in which the lick signal have best separation across conditions and will be analyzed
lickTimeRelativeEvent = 'reward';           
monkeyName = 'Mojo';  % 'Moj' or 'Mac'
groupingCriteria       = 'type' ; % this parameter shows based on what criteria trials are grouped.
Fs = 1000;   % assume 1K sampling rate
whichFigUWannaPlot = '5';  % '5' or '6'

%% Read out list of all the files related to this session
dataPath = [main_folder data_folder original_data_folder];
allFiles = dir(dataPath);    

% remove the '.' and '..' from the folders list and filter sessions for monkey names
idx = [];
for ses = 1:numel(allFiles) % for each session
    if (numel(allFiles(ses).name) < 3) || ~(strcmp(allFiles(ses).name(1:3), monkeyName(1:3)))  
        idx = [idx; ses];
    end
end
allFiles(idx) = [];
clear ses idx

%% Extract and group licking signal and rea
allSesGroupedData = {};

for ses = 1:numel(allFiles) % for each session
    cfg = [];
    sessionFolder = allFiles(ses).name;
 if 0 % ---------------------------------------------------------------------------------------------------------------------   
    bhvFileInfo = dir([dataPath sessionFolder '\*.bhv']); 
    bhvFilePath = [dataPath sessionFolder '\' bhvFileInfo.name];
        
    lfpFileInfo    = dir([dataPath sessionFolder '\*.ns2']);   
    if length(lfpFileInfo) > 1
        warning(['Multiple NS files in folder ' [dataPath sessionFolder]]);
        continue;
    end       
    lfpFilePath = [dataPath sessionFolder '\' lfpFileInfo.name];

    % form the configuration structure for 'ft_definetrial()'
    tmp_cfg                 = [];
    tmp_cfg.headerformat    = 'blackrock_nsx';
    tmp_cfg.data_dir        = dataPath;
    tmp_cfg.trialfun        = 'fieldtrip_trialfun_RiskBhv'; % string with function name, see below (default = 'ft_trialfun_general')
    tmp_cfg.trialdef.pre    = 1; % time in seconds; this only indicates absolut value, do not use '-'
    tmp_cfg.trialdef.post   = 1; % time in seconds
        
    tmp_cfg.headerfile  = lfpFilePath;
    tmp_cfg.dataset     = lfpFilePath;
    tmp_cfg.session_dir = sessionFolder;
        
    % Config for aligning to cue
    tmp_cfg.trialdef.eventtype  = alignEvent;
    tmp_cfg.trialdef.interval   = evetTimeInterval;

    % read out the data the same way it is done in ReadRiskDataBhv.m
    new_cfg = ft_definetrial(tmp_cfg);
    
    % save the extracted data on the hard drive for future use
    save([save_dir allFiles(ses).name '(' alignEvent ')'], 'new_cfg', 'alignEvent', 'evetTimeInterval')
 end % ----------------------------------------------------------------------------------------------------------------------

    load ([save_dir 'Bhv_' allFiles(ses).name])
 
    % group trials for 9 task conditions
    [output] = GroupTrials(new_cfg.event, groupingCriteria);

    % convert the event tructure to table
    eventTable = struct2table(new_cfg.event);

    % get the trial index for successful trials per condition
    for gr = 1:length(output)
        output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
        % for successful trials of this group, collect the lick signal in the asked time interval 
        tempLick = [];
        tempRT = [];
        bestSampleTimeStamp = []; % this is the time stamp in which the lick signal will be analyzed
        for tr = 1:length(output(gr).SuccessTrlIdx)
            trIdx = output(gr).SuccessTrlIdx(tr);
            if (ses==9 && trIdx==133 && strcmp(monkeyName(1:3),'Mac'))
                continue
            end
            % read out the lick signal at the predefined time point
            switch lickTimeRelativeEvent
                case 'cue'
                    bestSampleTimeStamp(tr,:) = eventTable.ActualEventTime(trIdx).CueOnset + (lickBestTime*Fs);
                    tmpSig = (eventTable.LickSignal{trIdx} > 2.5)';
                    tempLick = [tempLick; mean(tmpSig(bestSampleTimeStamp(tr,1):bestSampleTimeStamp(tr,2)))];
                    clear tmpSig
                case 'reward'  
                    bestSampleTimeStamp(tr,:) = eventTable.RewardOnTime{trIdx}(1) + (lickBestTime*Fs);
                    tmpSig = (eventTable.LickSignal{trIdx} > 2.5)';
                    tempLick = [tempLick; mean(tmpSig(bestSampleTimeStamp(tr,1):bestSampleTimeStamp(tr,2)))];
                    clear tmpSig
            end
            % read out the reaction time for this successful trial of this group
            tempRT = [tempRT; eventTable.ReactionTime(trIdx)];
            clear trIdx
        end
        output(gr).SuccessLick = tempLick; 
        output(gr).OriginalCutSampleLick = bestSampleTimeStamp;
        output(gr).SuccessRT   = tempRT; 
        clear tempLick tempRT tr bestSampleTimeStamp

    end
    
    allSesGroupedData{ses} = output;
    clear gr output eventTable new_cfg tmp_cfg temp 
end
clear ses

%%  plot either fig5 or fig6

eval(['plotFig' whichFigUWannaPlot])
