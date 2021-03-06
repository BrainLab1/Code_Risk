
% this script analyzes effect of the prior trial reward on the licking
% probability of the current trial

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

           
% % alignEvent = 'cue';           
eventTimeInterval = [-0.6, 0.1];
alignEvent = 'reward';           
% evetTimeInterval = [-0.55, 0.1];
monkeyName = 'Mojo';  % 'Moj' or 'Mac'
groupingCriteria       = 'PreTrlOutcome & CurrentTrlEV' ; % 'PreTrialEV' or 'PreTrialVAR'
Fs = 1000;   % assume 1K sampling rate

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
%%
allSesGroupedData = {};

for ses = 1:numel(allFiles) % for each session
    cfg = [];
    sessionFolder = allFiles(ses).name;

    load ([save_dir 'Bhv_' allFiles(ses).name])
 
    % group trials for 18 task conditions
    [output] = GroupTrials(new_cfg.event, groupingCriteria);   
    
    % if doesn't exist, make a new variable with the same size as the number of groups
    if ~exist('allSesSuccessLick','var')
        allSesSuccessLick = cell(length(output), 1);
    end

    % convert the event tructure to table
    eventTable = struct2table(new_cfg.event);
    
    % get the trial index for successful trials per condition
    for gr = 1:length(output)
        output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
        % for successful trials of this group, collect the lick signal in the asked time interval 
        temp = [];
        cutInterval = [];
        for tr = 1:length(output(gr).SuccessTrlIdx)
            trIdx = output(gr).SuccessTrlIdx(tr);
            if (ses==9 && trIdx==133)
                continue
            end
            switch alignEvent
                case 'cue'
                    cutInterval(tr,:) = eventTable.ActualEventTime(trIdx).CueOnset + (eventTimeInterval*Fs);
                    tmpSig = eventTable.LickSignal{trIdx} > 2.5;
                    temp = [temp; tmpSig(cutInterval(tr,1) : cutInterval(tr,2))'];
                    clear tmpSig
                case 'reward'
                    cutInterval(tr,:) = eventTable.RewardOnTime{trIdx}(1) + (eventTimeInterval*Fs);
                    tmpSig = eventTable.LickSignal{trIdx} > 2.5;
                    temp = [temp; tmpSig(cutInterval(tr,1) : cutInterval(tr,2))'];
                    clear tmpSig
            end
            clear trIdx
        end
        output(gr).SuccessLick = temp; 
        output(gr).OriginalCutSampleInterval = cutInterval;
        clear temp
%         temp = 
        allSesSuccessLick(gr) = {[allSesSuccessLick{gr}; output(gr).SuccessLick]};
        clear tr cutInterval
    end
    allSesGroupedData{ses} = output;
    clear gr output eventTable new_cfg tmp_cfg temp 
end
clear ses

%% plot each condition in its predefined position
xLimit = eventTimeInterval;
% yLimit = [0 550];
% plotBoxAspectRatio = [2, 1, 1];

figure, hold on, box on
title([monkeyName ', Alignment: ' alignEvent ', Conditions: ' regexprep(groupingCriteria, '_', '-')])

txtContent = {};
for gr = 1:length(allSesSuccessLick)
    varargout = shadedErrorBar( 1:size(allSesSuccessLick{gr},2),...
                                mean(allSesSuccessLick{gr}),...
                                std(allSesSuccessLick{gr})/sqrt(size(allSesSuccessLick{gr},1)),...
                                {'Color',colorCode(gr,:),'LineWidth',3},1);
          
    txtContent{gr} = ['\color[rgb]{' num2str(colorCode(gr,:)) '} ' num2str( getfield(allSesGroupedData{1}(gr), 'Value') )];      
end
text(mean(xlim), mean(ylim), txtContent)

switch alignEvent
    case 'cue'
        yl = [0 1]; %ylim;
        xTickPos = (-eventTimeInterval(1)*Fs) + [0 400 800];
        plot(repmat(xTickPos(1),1,2), yl, '-- k')
        plot(repmat(xTickPos(2),1,2), yl, '-- k')
        plot(repmat(xTickPos(3),1,2), yl, '-- k')
        xTickLabel = {'0', '0.4', '0.8'};
        text(xTickPos(1), 0.8*yl(2), 'CueOn', 'Rotation', 90)
        text(xTickPos(2), 0.8*yl(2), 'CueOff', 'Rotation', 90)
        text(xTickPos(3), 0.8*yl(2), 'TargetOn', 'Rotation', 90)
        set(gca, 'XTick', xTickPos, 'XTickLabel', xTickLabel)
        
    case 'reward'
        yl = [0 1];
        xTickPos = (-eventTimeInterval(1)*Fs) + [-377 0];
        plot(repmat(xTickPos(1),1,2), yl, '-- k')
        plot(repmat(xTickPos(2),1,2), yl, '-- k')
        xTickLabel = {'-0.377', '0'};
        text(xTickPos(1), 0.8*yl(2), 'TargetAcquired', 'Rotation', 90)
        text(xTickPos(2), 0.8*yl(2), 'RewardOn', 'Rotation', 90)
        set(gca, 'XTick', xTickPos, 'XTickLabel', xTickLabel)
        
end

xlabel(['Time(sec) aligned to ' alignEvent])
ylabel('Averaged licking signal')
