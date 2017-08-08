
clear all;
close all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
main_folder          = 'Z:\'; 
data_folder          = 'data\'; % 
original_data_folder = 'Risk\Original Data_Extracted\';
save_dir             = 'X:\Report 2\Bahareh\BehaviorData\';

%% Read out list of all the files related to this session
dataPath = [main_folder data_folder original_data_folder];
allFiles = dir(dataPath);    

% remove the '.' and '..' from the folders list and filter sessions for monkey names
idx = [];
for ses = 1:numel(allFiles) % for each session
    if numel(allFiles(ses).name) < 3  
        idx = [idx; ses];
    end
end
allFiles(idx) = [];
clear ses idx

%% Extract and group licking signal and read
% 
allSesEvtTime = struct('MLTrialStart', [], 'TrialGateOn', [], 'FixationOn', [], 'FixationAcquired', [], ...
                        'CueOnset', [], 'CueOffset', [], 'TargetOnset', [], 'TargetAcquired', [], ...
                        'TargetOffset', [], 'TrialGateOff', [], 'MLTrialEnd', [], 'FirstRewardPulse', [], 'RewardEndTime', []);

% get the list of the event names in correct trials
evtNames = fieldnames(allSesEvtTime);
ind = 1;
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

    load ([save_dir allFiles(ses).name '(cue)'])
    
    % convert the structure arry to a table
    eventTable = struct2table(new_cfg.event);
    
    % find unsuccessful trials for this session
    trlIdx = find(eventTable.TrialErrorCode ~= 0);
    
    % remove the error trials from the main event stucture
    new_cfg.event(trlIdx) = [];
    clear eventTable
    
    % convert the structure arry to a table but this time it only includes correct trials
    eventTable = struct2table(new_cfg.event);
    
    % delet the new_cfg structure
    clear new_cfg
    
    % 
    c = struct2cell(eventTable.ActualEventTime)';
    
    % collect the event times for this session
    ev = 1;
    while ev <= length(evtNames)  
        switch ev
            case 1 % this is for MLTrialStart
                aa = [];
                for tr=1:size(eventTable,1)
                    aa = [aa; mean(eventTable.ActualEventTime(tr).MLTrialStart)];
                end      
                allSesEvtTime.(evtNames{ev}) = [allSesEvtTime.(evtNames{ev}); aa];
                clear aa
                
            case 11 % this is for MLTrialEnd
                aa = [];
                for tr=1:size(eventTable,1)
                    aa = [aa; mean(eventTable.ActualEventTime(tr).MLTrialEnd)];
                end      
                allSesEvtTime.(evtNames{ev}) = [allSesEvtTime.(evtNames{ev}); aa];
                clear aa
                
            case 12 % this is for FirstRewardPuls
                aa = [];
                for tr=1:size(eventTable,1)
                    aa = [aa; eventTable.RewardOnTime{tr}(1)];
                end      
                allSesEvtTime.(evtNames{ev}) = [allSesEvtTime.(evtNames{ev}); aa];
                clear aa
                
            case 13 % this is for RewardEndTime
                aa = [];
                for tr=1:size(eventTable,1)
                    aa = [aa; eventTable.RewardOffTime{tr}(end)];
                end      
                allSesEvtTime.(evtNames{ev}) = [allSesEvtTime.(evtNames{ev}); aa];
                clear aa

            otherwise
                allSesEvtTime.(evtNames{ev}) = [allSesEvtTime.(evtNames{ev}); cell2mat(c(:,ev))];
        end
        ev = ev+1;
    end
    
    clear eventTable
    
    if isempty(strfind(sessionFolder,'Mojo'))
        allSesEvtTime.Monkey(ind:ind+tr,1) = 1;
    else
        allSesEvtTime.Monkey(ind:ind+tr,1) = 2;
    end
    ind = tr + ind;
    
end
clear ses

%% plot the distribution of inter-event time
numRows = 4;
numColumns = 3;
binWidth = 5; % bin width in ms

monkey_ind1 = find(allSesEvtTime.Monkey == 1);
monkey_ind2 = find(allSesEvtTime.Monkey == 2);
monkey_ind2(end) = [];

figure('Name', 'Inter-event timing')
ev = 2;
while ev <= length(evtNames)-1  % exclude TrialGateOn and TrialGateOff
    subplot(numRows, numColumns, ev-1), hold on, box on
    title([evtNames{ev} ' - ' evtNames{ev-1}])
    xlabel('Time in (ms)')
    ylabel('num. trials')
    
    % for Mojo
    tmp1 = allSesEvtTime.(evtNames{ev});
    tmp2 = allSesEvtTime.(evtNames{ev-1});
    deltaTime = tmp1(monkey_ind2) - tmp2(monkey_ind2);
    xLimit = [min(deltaTime)  max(deltaTime)] + [-10  10];
    xbins = xLimit(1) : binWidth : xLimit(2);
    [counts,~] = hist(deltaTime, xbins);
    bar(xbins, counts, 'FaceColor', [153 255 204]/255)
    set(gca, 'XLim', xLimit)
    clear tmp1 tmp2
    
    % for MacDuff
    tmp1 = allSesEvtTime.(evtNames{ev});
    tmp2 = allSesEvtTime.(evtNames{ev-1});
    deltaTime = tmp1(monkey_ind1) - tmp2(monkey_ind1);
    xLimit = [min(deltaTime)  max(deltaTime)] + [-10  10];
    xbins = xLimit(1) : binWidth : xLimit(2);
    [counts,~] = hist(deltaTime, xbins);
    bar(xbins, counts, 'FaceColor', [153 100 204]/255)
    set(gca, 'XLim', xLimit)
    clear tmp1 tmp2
    
    ev = ev+1;
end

    subplot(numRows, numColumns, 11), hold on, box on
    title([evtNames{12} ' - ' evtNames{8}])
    xlabel('Time in (ms)')
    ylabel('num. trials')
    
    % for Mojo
    tmp1 = allSesEvtTime.(evtNames{12});
    tmp2 = allSesEvtTime.(evtNames{8});
    deltaTime = tmp1(monkey_ind2) - tmp2(monkey_ind2);
    xLimit = [min(deltaTime)  max(deltaTime)] + [-10  10];
    xbins = xLimit(1) : binWidth : xLimit(2);
    [counts,~] = hist(deltaTime, xbins);
    bar(xbins, counts, 'FaceColor', [153 255 204]/255)
    set(gca, 'XLim', xLimit)
    clear tmp1 tmp2

    % for MacDuff
    tmp1 = allSesEvtTime.(evtNames{12});
    tmp2 = allSesEvtTime.(evtNames{8});
    deltaTime = tmp1(monkey_ind1) - tmp2(monkey_ind1);
    xLimit = [min(deltaTime)  max(deltaTime)] + [-10  10];
    xbins = xLimit(1) : binWidth : xLimit(2);
    [counts,~] = hist(deltaTime, xbins);
    bar(xbins, counts, 'FaceColor', [153 100 204]/255)
    set(gca, 'XLim', xLimit)
    clear tmp1 tmp2
    
    subplot(numRows, numColumns, 12), hold on, box on
    title([evtNames{10} ' - ' evtNames{13}])
    xlabel('Time in (ms)')
    ylabel('num. trials')
    
    % for Mojo
    tmp1 = allSesEvtTime.(evtNames{10});
    tmp2 = allSesEvtTime.(evtNames{13});
    deltaTime = tmp1(monkey_ind2) - tmp2(monkey_ind2);
    xLimit = [min(deltaTime)  max(deltaTime)] + [-10  10];
    xbins = xLimit(1) : binWidth : xLimit(2);
    [counts,~] = hist(deltaTime, xbins);
    bar(xbins, counts, 'FaceColor', [153 255 204]/255)
    set(gca, 'XLim', xLimit)
    clear tmp1 tmp2
    
    % for MacDuff
    tmp1 = allSesEvtTime.(evtNames{10});
    tmp2 = allSesEvtTime.(evtNames{13});
    deltaTime = tmp1(monkey_ind1) - tmp2(monkey_ind1);
    xLimit = [min(deltaTime)  max(deltaTime)] + [-10  10];
    xbins = xLimit(1) : binWidth : xLimit(2);
    [counts,~] = hist(deltaTime, xbins);
    bar(xbins, counts, 'FaceColor', [153 100 204]/255)
    set(gca, 'XLim', xLimit)
    
