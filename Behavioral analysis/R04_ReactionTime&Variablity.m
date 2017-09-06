
% this script plots averaged reaction time by current trial expected reward
% and variance; for data readout and analysis followings considered:
% - remove outlier trials based on outlier indices saved with the event structure  
% - use DiodeReactionTime (i.e., DRT)
% - transform DRT to zscore(log(DRT))
% - run non-parametric anova using ART

clear all;
clc;
dbstop if error
format short

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
data_folder             = 'Z:\data\Risk\Behavior\';
csv_save_path           = 'D:\Report04_ART for RT\';

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
groupingCriteria       = 'expecter_reward & RewardVariance' ; % this parameter shows based on what criteria trials are grouped.
Fs = 1000;   % assume 1K sampling rate

%% add path
addpath(genpath('D:\Code_localGit'))

%% Read out list of all the files related to this session
dataPath = [data_folder];
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
%% Extract, group and transform diode reaction time
allSesGroupedData = {};

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
    
    thisSesDRT = [];
    % get the trial index for successful trials per condition
    for gr = 1:length(output)
        output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
        numSuccInGrp(gr,1) = length(output(gr).SuccessTrlIdx);
        tempDRT = [];
        for tr = 1:length(output(gr).SuccessTrlIdx) % for each success trial of this group
            trIdx = output(gr).SuccessTrlIdx(tr);
            % read out the diode-reaction time for this successful trial of this group
            tempDRT = [tempDRT; eventTable.DiodeReactionTime(trIdx)];
            clear trIdx
        end
        output(gr).SuccessDRT   = tempDRT; 
        clear tempDRT tr 
        thisSesDRT = [thisSesDRT; output(gr).SuccessDRT];
    end
    clear gr 
    
    % if there are 0 RTs, add eps to them to prevent -Inf values after  log() function
    thisSesDRT(thisSesDRT == 0) = eps;
    % transform diode reaction times for this session
    trnsfrmDRT = zscore(log(thisSesDRT));
    aa = mat2cell(trnsfrmDRT, numSuccInGrp, 1);
    % add the transformed DRT for each group
    for gr = 1:length(output)
        output(gr).TrnsfrmDRT = aa{gr};
        output(gr).TrnsfrmDRTvar = var(aa{gr});
    end
    clear aa gr trnsfrmDRT thisSesDRT
    
    allSesGroupedData{ses} = output;
    clear output eventTable new_cfg numSuccInGrp
end
clear ses

%% collect the grouped DRTs from all sessions 
allSesTrnsfrmDRT = []; %collected RT and condition info from all sessions -> [TrnsfrmDRT, expReward, rewardVar sessionID]
allSesDRTVar     = []; %collected RT variance for each session and condition -> [Var(TrnsfrmDRT), expReward, rewardVar sessionID]
for ses = 1:length(allSesGroupedData)
    expReward = [];
    rewardVar = [];
    for gr=1:length(allSesGroupedData{ses})
        expReward(gr) = allSesGroupedData{ses}(gr).Value(1);
        rewardVar(gr) = allSesGroupedData{ses}(gr).Value(2);
        allSesTrnsfrmDRT = [ allSesTrnsfrmDRT; ...
                             [allSesGroupedData{ses}(gr).TrnsfrmDRT, ...
                              repmat(expReward(gr), length(allSesGroupedData{ses}(gr).TrnsfrmDRT), 1), ...
                              repmat(rewardVar(gr), length(allSesGroupedData{ses}(gr).TrnsfrmDRT), 1)], ...
                              repmat(ses, length(allSesGroupedData{ses}(gr).TrnsfrmDRT), 1)];
       
        % collect RT variablity for this session                  
        allSesDRTVar = [ allSesDRTVar; ...
                         [var(allSesGroupedData{ses}(gr).TrnsfrmDRT), expReward(gr), rewardVar(gr), ses] ];                           
    end
    clear gr
end
clear ses 

%% Apply ART for reaction times
% write the csv file and start from the second row
csvFileName = [csv_save_path monkeyName '_RT.csv']; % full path plus file name for writing data in the csv file
csvVarFileName = [csv_save_path monkeyName '_RTvar.csv']; % full path plus file name for writing RT variance data in the csv file
% I had to round the close-to-zero values otherwise ART was returning error
aa = allSesTrnsfrmDRT(:,1);
aa(find((-0.001<aa) .* (aa<0.001))) = 0;
csvwrite(csvFileName, [ [-1000 -1111 -2222 -3333 ]; [ones(size(allSesTrnsfrmDRT,1),1), allSesTrnsfrmDRT(:,2:3), aa] ], 0, 0);
clear aa
% write RT variance data to the csv file
csvwrite(csvVarFileName, [ [-1000 -1111 -2222 -3333 ]; [ones(size(allSesDRTVar,1),1), allSesDRTVar(:,2:3), allSesDRTVar(:,1)] ], 0, 0);

% ......... Apply ART manually and get the .rt.csv file .........
display('****  Run ART manually!, then press any key to continue  ****')
pause
% readout the aligned-tanked RT data
artData = csvread([csvFileName(1:end-4) '.art.csv'],1); % .art data -> [subjID EV Var Data, A(EV) A(Var) A(interaction), AR(EV) AR(Var) AR(interaction)] (row/column indices start from 0)
txtRT = [{'non-parametric ANOVA'};{'--------------------'}];
% run parametric anova on A.R. data for EV
[~,tbl_EV] = anovan(artData(:,8), artData(:,2:3),'model','interaction','varnames',{'EV','Var'},'display','off');
txtRT = [txtRT; {['p-' tbl_EV{2,1} '= ' num2str(tbl_EV{2,7})]}];
% run parametric anova on A.R. data for Var
[~,tbl_Var] = anovan(artData(:,9), artData(:,2:3),'model','interaction','varnames',{'EV','Var'},'display','off');
txtRT = [txtRT; {['p-' tbl_Var{3,1} '= ' num2str(tbl_Var{3,7})]}];
% run parametric anova on A.R. data for Interaction
[~,tbl_Int] = anovan(artData(:,10), artData(:,2:3),'model','interaction','varnames',{'EV','Var'},'display','off');
txtRT = [txtRT; {['p-' tbl_Int{4,1} '= ' num2str(tbl_EV{4,7})]}];
clear tbl_EV tbl_Var tbl_Int

% readout the aligned-tanked RT data
artData = csvread([csvVarFileName(1:end-4) '.art.csv'],1); % .art data -> [subjID EV Var Data, A(EV) A(Var) A(interaction), AR(EV) AR(Var) AR(interaction)] (row/column indices start from 0)
txtRTvar = [{'non-parametric ANOVA'};{'--------------------'}];
% run parametric anova on A.R. data for EV
[~,tbl_EV] = anovan(artData(:,8), artData(:,2:3),'model','interaction','varnames',{'EV','Var'},'display','off');
txtRTvar = [txtRTvar; {['p-' tbl_EV{2,1} '= ' num2str(tbl_EV{2,7})]}];
% run parametric anova on A.R. data for Var
[~,tbl_Var] = anovan(artData(:,9), artData(:,2:3),'model','interaction','varnames',{'EV','Var'},'display','off');
txtRTvar = [txtRTvar; {['p-' tbl_Var{3,1} '= ' num2str(tbl_Var{3,7})]}];
% run parametric anova on A.R. data for Interaction
[~,tbl_Int] = anovan(artData(:,10), artData(:,2:3),'model','interaction','varnames',{'EV','Var'},'display','off');
txtRTvar = [txtRTvar; {['p-' tbl_Int{4,1} '= ' num2str(tbl_EV{4,7})]}];
clear tbl_EV tbl_Var tbl_Int
%% plot reaction times
% find trials with [ExpReward Variance] 
expReward = [3 6 9];
rewVar = [0 1 4];
legTxt = {};
fg1 = figure('Name', monkeyName);
hold on, box on
title(monkeyName)
xlabel('Expected reward')
ylabel({'Ave. Z(log( reaction time ))','+/- SEM'})

for vr=1:3
    legTxt = [legTxt; {['\bf \color[rgb]{' num2str(colorCode(vr,:)) '} Var= ' num2str(rewVar(vr))]}];
    for ex = 1:3
        idx = find( (allSesTrnsfrmDRT(:,2) == expReward(ex)) .* (allSesTrnsfrmDRT(:,3)==rewVar(vr)) );
        eval(['crv' num2str(rewVar(vr)) '(ex) = mean(allSesTrnsfrmDRT(idx,1));'])
        eval(['err' num2str(rewVar(vr)) '(ex) = std(allSesTrnsfrmDRT(idx,1))/sqrt(length(idx));'])
        clear idx
    end
    eval(['errorbar( [1:3], crv' num2str(rewVar(vr)) ', err' num2str(rewVar(vr)) ', ''LineWidth'', 2, ''Color'', colorCode(vr,:) )'])
end
clear vr ex

set(gca, 'XTick', [1:3], 'XTickLabel', {'3','6','9'}, 'XLim', [0 4])

text(0.2, 0.7*max(ylim), txtRT)
text(3.2, min(ylim)+0.05, legTxt)
clear legTxt crv0 err0 crv1 err1 crv4 err4

%% plot reactio time variance
% find trials with [ExpReward Variance] 
expReward = [3 6 9];
rewVar = [0 1 4];
legTxt = {};
fg2 = figure('Name', monkeyName);
hold on, box on
title(monkeyName)
xlabel('Expected reward')
ylabel({'Variance of Z(log( reaction time ))','+/- SEM'})

for vr=1:3
    legTxt = [legTxt; {['\bf \color[rgb]{' num2str(colorCode(vr,:)) '} Var= ' num2str(rewVar(vr))]}];
    for ex = 1:3
        idx = find( (allSesDRTVar(:,2) == expReward(ex)) .* (allSesDRTVar(:,3)==rewVar(vr)) );
        eval(['crv' num2str(rewVar(vr)) '(ex) = mean(allSesDRTVar(idx,1));'])
        eval(['err' num2str(rewVar(vr)) '(ex) = std(allSesDRTVar(idx,1))/sqrt(length(idx));'])
        clear idx
    end
    eval(['errorbar( [1:3], crv' num2str(rewVar(vr)) ', err' num2str(rewVar(vr)) ', ''LineWidth'', 2, ''Color'', colorCode(vr,:) )'])
end
clear vr ex

set(gca, 'XTick', [1:3], 'XTickLabel', {'3','6','9'}, 'XLim', [0 4])

text(0.2, 0.7*max(ylim), txtRTvar)
text(3.2, min(ylim)+0.05, legTxt)
