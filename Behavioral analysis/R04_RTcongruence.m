

% This script takes the reaction times;
% transforms them to zscore( Log(reaction time) ); removes ooutliers;
% analyzes effect of expected reward and congruence

clear all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
data_folder             = 'Z:\data\Risk\Behavior\';

colorCode = [  204 204 255  % F3
               255 153 153  % F6
               153 255 204   % F9
               102 102 255   % B3L
               255 51  51    % B6L
               51  255 102   % B9L              
               0   0   255   % B3H
               204 0   0     % B6H
               0  153 153]/255; % B9H   
           
monkeyName = 'Mojo';  % 'Moj' or 'Mac'
groupingCriteria       = 'expected_reward & Congruence' ; % this parameter shows based on what criteria trials are grouped.

%% add path
addpath(genpath('D:\Code_localGit'))

%% Read out list of all the files related to this session
dataPath = data_folder;
allFiles = dir(dataPath);    

% remove the '.' and '..' from the folders list and filter sessions for monkey names
idx = [];
for ses = 1:numel(allFiles) % for each session
    if (numel(allFiles(ses).name) < 3) || ~(strcmp(allFiles(ses).name(5:7), monkeyName(1:3)))  
        idx = [idx; ses];
    end
end
allFiles(idx) = [];
clear ses idx

%% Extract and group reaction times
allSesGroupedData = {};
allSesSuccRT = cell(numel(allFiles),1);  % collected reaction times from successful trials of each session

for ses = 1:numel(allFiles) % for each session
    load ([dataPath allFiles(ses).name])
    
    % convert the event tructure to table
    eventTable = struct2table(new_cfg.event);
    
    % ........remove outliers....
    % find outlier trials for duration of cue presentation & reward time and reaction time
    outTrials = find(eventTable.OutCueOnOffZ2 + eventTable.OutAqToRwdZ2 + eventTable.OutRTZ2);
    clear eventTable
    % remove oulier trials from eventTable and new_cfg.event 
    new_cfg.event(outTrials) = [];
    clear outTrials
 
    % group trials for 27 task conditions
    [output] = GroupTrials(new_cfg.event, groupingCriteria);   

    % convert the outlier-removed event tructure to table
    eventTable = struct2table(new_cfg.event);
    
    % get the trial index for successful trials per condition
    for gr = 1:length(output)
        output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
        % for successful trials of this group, collect the RT 
        tempRT = [];
        for tr = 1:length(output(gr).SuccessTrlIdx)
            trIdx = output(gr).SuccessTrlIdx(tr);
            % read out the reaction time for this successful trial of this group
            tempRT = [tempRT; eventTable.DiodeReactionTime(trIdx)];
            clear trIdx
        end
        % if there are 0 RTs, add eps to them to prevent -Inf values after  log() function
        tempRT(tempRT == 0) = eps;
        % add diode reaction times to the output
        output(gr).SuccessRT   = tempRT;
        allSesSuccRT(ses) = {[allSesSuccRT{ses}; output(gr).SuccessRT]};
        clear  tempRT tr 
    end
        
    

    allSesGroupedData{ses} = output;
    clear gr output eventTable new_cfg tmp_cfg temp 
end
clear ses

%% take the zscore(log(RT)), remove outliers
allSesSuccTransformRT = cell(numel(allFiles),1);  % transformed reaction times from successful trials of each session
for ses = 1:length(allSesSuccRT)
    % transform reaction times
    allSesSuccTransformRT{ses} = zscore(log(allSesSuccRT{ses}));
    
    trlPerGrp = []; % number of successful trials per group for this session
    tmp = [];
    for gr = 1:length(allSesGroupedData{ses})
        trlPerGrp(gr,1) = length(allSesGroupedData{ses}(gr).SuccessRT);
        tmp = [tmp; allSesGroupedData{ses}(gr).SuccessRT]; % collect reaction times in order of groups
    end
    clear gr
    % transform the collected reaction times for this session
    tmp = zscore(log(tmp));
    tmp = mat2cell(tmp, trlPerGrp, 1);
    % add the transformed RTs for to each group
    for gr = 1:length(allSesGroupedData{ses})
        allSesGroupedData{ses}(gr).TransformSuccRT = tmp{gr};
    end
    clear tmp
end

%% plot distribution of reaction times from all sessions
figure('Name', monkeyName)
hold on, box on, 
title([monkeyName ', ' num2str(length(allSesGroupedData)) ' sessions, ' ...
       num2str(length(cell2mat(allSesSuccTransformRT))) ' trials'])
ylabel('number of trials')                    
xlabel('zscore( Log( ReactionTime ) )')
aa = cell2mat(allSesSuccTransformRT);  
xLimit = [min(aa) max(aa)];
xbin = xLimit(1)-0.5 : 0.1 : xLimit(2)+0.5;
counts = hist(aa , xbin);
bar(xbin, counts, 'FaceColor', colorCode(3,:))
yLimit = ylim;
plot([-2 -2], ylim, 'r--')
plot([2 2], ylim, 'r--')
plot([0 0], ylim, 'k--')
set( gca, 'XTick', sort(unique([get(gca, 'XTick') -2 2])) )
clear aa counts

%% RT (cong-ic trials), separate by EV (but collapse across variances) and use a 1-way KW analysis 
meanRT = [];  % groups ->  [[3 1], [6 1], [9 1], [3 -1], [6 -1], [9 -1]]
for ses = 1:length(allSesSuccRT)
    for gr = 1:length(allSesGroupedData{ses})
        meanRT(ses, gr) = mean(allSesGroupedData{ses}(gr).TransformSuccRT);
    end
end
clear ses gr

% calculate the congRT-incongRT
deltaRT = meanRT(:,1:3) - meanRT(:,4:6);        

%% plot the result
figure, hold on, box on
xlabel('expected reward')
ylabel([{'congRT - incongRT'},{'+/- SEM'}])
title([monkeyName ', ' num2str(length(allSesSuccRT)) ' sessions'])
errorbar(1:3, mean(deltaRT), std(deltaRT)/sqrt(size(deltaRT,1)), 'LineWidth', 2)
set(gca, 'XLim', [0 4], 'XTick', [1:3], 'XTickLabel', [{'3'},{'6'},{'9'}])
[pKruskalWallis,tbl,stats] = kruskalwallis(deltaRT,[],'off');
text(mean(xlim), mean(ylim), ['pKruskalWallis = ' num2str(pKruskalWallis)])
                    
                    
                    
                    
                    
