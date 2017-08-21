
% this script analyzes effect of the prior trial reward on the licking
% probability of the current trial

clear
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
monkeyName = 'MacDuff';  % 'Moj' or 'Mac'
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
indx = 1;

for ses = 1:numel(allFiles) % for each session
    cfg = [];
    sessionFolder = allFiles(ses).name;

    load ([save_dir 'Bhv_' allFiles(ses).name])
    
    eventTable = new_cfg.event;
    for tr = 2:size(eventTable,1)  % for each trial
        if ~eventTable(tr).TrialErrorCode % if current trial was successful
            if ( ~eventTable(tr-1).TrialErrorCode  && ~strcmp(eventTable(tr-1).type,'F') )
                temp = [];
                switch alignEvent
                    case 'cue'
                        cutInterval(tr,:) = eventTable.ActualEventTime(tr).CueOnset + (eventTimeInterval*Fs);
                        tmpSig = eventTable.LickSignal{tr} > 2.5;
                        temp = [tmpSig(cutInterval(tr,1) : cutInterval(tr,2))'];
                        clear tmpSig
                    case 'reward'
                        cutInterval(tr,:) = eventTable(tr).RewardOnTime(1) + (eventTimeInterval*Fs);
                        tmpSig = eventTable(tr).LickSignal > 2.5;
                        temp = [tmpSig(cutInterval(tr,1) : cutInterval(tr,2))'];
                        clear tmpSig
                end
                y(indx) = mean(temp);
                
                X(indx,1) = eventTable(tr).expected_reward;
                X(indx,2) = eventTable(tr).RewardVariance;
                X(indx,3) = eventTable(tr-1).RewardVariance;
                X(indx,4) = eventTable(tr-1).TotalRewardTime;
                
                if eventTable(tr-1).TotalRewardTime < ((eventTable(tr-1).expected_reward/3-1)*100+125)% low rewarded pre-trial
                    X(indx,5) = -1;
                else % high rewarded pre-trial
                    X(indx,5) = 1;
                end
                indx = indx + 1;
            end
        end
    end
     
end
clear ses


[b,dev,stats] = glmfit(X,y,'normal');
    


