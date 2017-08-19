

% created on 12.08.2017 by Bahareh: This script takes the reaction times;
% transforms them to zscore( Log(reaction time) ); removes ooutliers;
% analyzes effect of expected reward and congruence

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
removeNegativeRT = 1; % this is a flag
monkeyName = 'Mac';  % 'Moj' or 'Mac'
groupingCriteria       = 'expected_reward & Congruence' ; % this parameter shows based on what criteria trials are grouped.

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

%% Extract and group reaction times
allSesGroupedData = {};
allSesSuccRT = cell(numel(allFiles),1);  % collected reaction times from successful trials of each session
removedTrials = {};

for ses = 1:numel(allFiles) % for each session
    sessionFolder = allFiles(ses).name;
    load ([save_dir 'Bhv_' allFiles(ses).name])
    
    % convert the event tructure to table
    eventTable = struct2table(new_cfg.event);
    
    % reamove trials with negative reaction times, if asked
    if removeNegativeRT
        invalidTrls = find(eventTable.DiodeReactionTime<=0);
        refineevent = new_cfg.event;
        refineevent(invalidTrls)=[];
        eventTable = struct2table(refineevent);
        removedTrials{ses} = invalidTrls;
        clear refineevent invalidTrls
    end

    % group trials for 9 task conditions
    [output] = GroupTrials(table2struct(eventTable), groupingCriteria);

    % get the trial index for successful trials per condition
    for gr = 1:length(output)
        output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
        % for successful trials of this group, collect the RT 
        tempRT = [];
        for tr = 1:length(output(gr).SuccessTrlIdx)
            trIdx = output(gr).SuccessTrlIdx(tr);
% % %             if (strcmp(monkeyName,'MacDuff') && ses==9 && trIdx==133)
% % %                 continue
% % %             end
            % read out the reaction time for this successful trial of this group
            tempRT = [tempRT; eventTable.DiodeReactionTime(trIdx)];
            clear trIdx
        end
        output(gr).SuccessRT   = tempRT; 
        allSesSuccRT(ses) = {[allSesSuccRT{ses}; output(gr).SuccessRT]};
        clear  tempRT tr 

    end
    
    allSesGroupedData{ses} = output;
    clear gr output eventTable new_cfg tmp_cfg temp 
end
clear ses

%% remove negative reaction times, take the zscore(log(RT)), remove outliers
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
if removeNegativeRT
    title([{[monkeyName ', ' num2str(length(allSesSuccRT)) ' sessions, ' num2str(length(cell2mat(allSesSuccTransformRT))) ' trials']},...
           {[num2str(length(cell2mat(removedTrials'))) ' trials with RT<=0 removed']}])
else
    title([monkeyName ', ' num2str(length(allSesGroupedData)) ' sessions, ' ...
           num2str(length(cell2mat(allSesSuccTransformRT))) ' trials'])
end
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

%% remove outliers from groups
for ses = 1:length(allSesSuccRT)
    for gr = 1:length(allSesGroupedData{ses})
        temp = allSesGroupedData{ses}(gr).TransformSuccRT;
        idx = find( (-rtOutlierFactor < temp) .* (temp < rtOutlierFactor) );
        allSesGroupedData{ses}(gr).TransformSuccRT = temp(idx);
        clear idx temp
    end
end
clear gr ses

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
if removeNegativeRT
    title([{[monkeyName ', ' num2str(length(allSesSuccRT)) ' sessions, outliers: std*' num2str(rtOutlierFactor)]},...
           {[num2str(length(cell2mat(removedTrials'))) ' trials with RT<=0 removed']}])
else
    title([monkeyName ', ' num2str(length(allSesSuccRT)) ' sessions, outliers: std*' num2str(rtOutlierFactor)])
end
errorbar(1:3, mean(deltaRT), std(deltaRT)/sqrt(size(deltaRT,1)))
set(gca, 'XLim', [0 4], 'XTick', [1:3], 'XTickLabel', [{'3'},{'6'},{'9'}])
[pKruskalWallis,tbl,stats] = kruskalwallis(deltaRT,[],'off');
text(mean(xlim), mean(ylim), ['pKruskalWallis = ' num2str(pKruskalWallis)])
                    
%% two factor analysis of variance with expected reward and congruence as factors                    
                    
allGroup = [];
for ses = 1:length(allSesSuccRT)
    for gr = 1:length(allSesGroupedData{ses})
        numTrl = length(allSesGroupedData{ses}(gr).TransformSuccRT);
        allGroup = [allGroup; [allSesGroupedData{ses}(gr).TransformSuccRT, repmat(allSesGroupedData{ses}(gr).Value, numTrl, 1)] ];
        clear numTrl
    end
end
clear ses gr

[p,tbl,stats] = anovan(allGroup(:,1), allGroup(:,2:3), 'model', 'interaction', 'varnames', {'ExpectedReward','Congruence'});


                    
                    
                    
                    
                    
                    
