
% this script plots and analyzes effect of the prior trial EV and Outcome on
% the licking probability of the current trial, in four time intervals 
% visual event times were considered based on the photo diode signal

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

           
eventTimeInterval = {[-0.3, 0]; [0 0.4]; [0.4 0.8]; [-0.377 0]};
alignEvent = {'cue'; 'cue'; 'cue'; 'reward'};           
monkeyName = 'Mojo';  % 'Moj' or 'Mac'
csv_save_path           = ['D:\Report04_ART for licking data\' monkeyName '\'];
groupingCriteria       = 'PreTrlEV & PreTrlOutcome & CurrTrlEV' ; 
Fs = 1000;   % assume 1K sampling rate

% each row of the 'allConditions' contains following values -> 
%                [PreTrlEV   PreTrlOutcome  CurrTrlEV] 
allConditions = [3  -1  3   %  -
                 3  -1  6   % | 
                 3  -1  9   %  -
                 3  0   3   %  -
                 3  0   6   % | 
                 3  0   9   %  -
                 3  1   3   %  -
                 3  1   6   % | 
                 3  1   9   %  -
                 ...
                 6  -1  3   %  -
                 6  -1  6   % | 
                 6  -1  9   %  -
                 6  0   3   %  -
                 6  0   6   % | 
                 6  0   9   %  -
                 6  1   3   %  -
                 6  1   6   % | 
                 6  1   9   %  -
                 ...
                 9  -1  3   %  -
                 9  -1  6   % | 
                 9  -1  9   %  -
                 9  0   3   %  -
                 9  0   6   % | 
                 9  0   9   %  -
                 9  1   3   %  -
                 9  1   6   % | 
                 9  1   9]; %  -                        
  
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
%%
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
    
    % if doesn't exist, make a new variable with the same size as the number of groups
    if ~exist('allSesSuccessLick','var')
        allSesSuccessLick = cell(length(output), length(alignEvent));
    end
    
    % convert the outlier-removed event tructure to table
    eventTable = struct2table(new_cfg.event);
    

    for tim = 1:length(alignEvent) % for each time interval
        for gr = 1:length(output)
            % get the trial index for successful trials per condition
            output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
            % get the total reward time for previous trials of output(gr).SuccessTrlIdx
            output(gr).PreTotalReward = eventTable.TotalRewardTime(output(gr).SuccessTrlIdx - 1);
            % get the reward variance for previous trials of output(gr).SuccessTrlIdx
            output(gr).PreVar = eventTable.RewardVariance(output(gr).SuccessTrlIdx - 1);
            % get the previous trial outcome 
            output(gr).PreOutcome = repmat(output(gr).Value(2), length(output(gr).TrialIdx),1);
            % get the previous trial expected value
            output(gr).PreEV = repmat(output(gr).Value(1), length(output(gr).TrialIdx),1);
            % get the reward variance for output(gr).SuccessTrlIdx trials
            output(gr).CurrVar = eventTable.RewardVariance(output(gr).SuccessTrlIdx);
            % get the expected reward for output(gr).SuccessTrlIdx trials
            output(gr).CurrEV = eventTable.expected_reward(output(gr).SuccessTrlIdx);

            % for successful trials of this group, collect the lick signal in the asked time interval 
            temp = [];
            cutInterval = [];
            for tr = 1:length(output(gr).SuccessTrlIdx)
                trIdx = output(gr).SuccessTrlIdx(tr);
                switch alignEvent{tim}
                    case 'cue'
                        % cut the licking signall based on DiodeEventTimes
                        cutInterval(tr,:) = eventTable.DiodeEventTime{trIdx}.CueOnset + (eventTimeInterval{tim}*Fs);
                        tmpSig = eventTable.LickSignal{trIdx} > 2.5;
                        temp = [temp; tmpSig(cutInterval(tr,1) : cutInterval(tr,2))'];
                        clear tmpSig
                    case 'reward'
                        cutInterval(tr,:) = eventTable.RewardOnTime{trIdx}(1) + (eventTimeInterval{tim}*Fs);
                        tmpSig = eventTable.LickSignal{trIdx} > 2.5;
                        temp = [temp; tmpSig(cutInterval(tr,1) : cutInterval(tr,2))'];
                        clear tmpSig
                end
                clear trIdx
            end
            eval(['output(gr).SuccessLickInterval' num2str(tim) ' = temp;'])
            eval(['output(gr).OrigCutSampleInterval' num2str(tim) ' = cutInterval;'])
            clear temp

            % collect lick signals of this group in this session's variable
            eval(['allSesSuccessLick(gr,tim) = {[allSesSuccessLick{gr,tim}; output(gr).SuccessLickInterval' num2str(tim) ']};'])
            clear tr cutInterval
        end
    end
    allSesGroupedData{ses} = output;
    clear gr output eventTable new_cfg tmp_cfg temp tim
end
clear ses

%% plot each condition in its predefined position
plotRowCurrEV = [3; 6; 9];
figure('Name', monkeyName)
numRows = 3; % each row is one current-trial-expected-value
numColumns = 4; % each column is one time interval
preEV = [3;6;9];
preOut = [-1;0;1];  % previouse-trila outcome [lose neutral win]
pltCounter = 0;

for rw = 1:numRows
    for clm = 1:numColumns
        pltCounter = pltCounter+1;
        fg = subplot(numRows, numColumns, pltCounter); hold on, box on
        txtContent = {};  % this is the text which will show the legend
        csvOutput = [];  % this is the matrix which will be written in the csv file
        csvFileName = [csv_save_path monkeyName num2str(rw) num2str(clm) '.csv']; % full path plus file name for writing data of this frame in the csv file
        for prot = 1:length(preOut)  %for each preOut
            currCrv = [];  % current curve!
            currErr = [];  % current error
            % form the curve with thrwee points 
            for prev = 1:length(preEV) % for each preEV
                % find conditions to be ploted in this row
                cdIdx = find( (allConditions(:,3) == plotRowCurrEV(rw)) .* (allConditions(:,2) == preOut(prot)) .* (allConditions(:,1) == preEV(prev)) );
                currCrv(1, prev) = mean(mean(allSesSuccessLick{cdIdx,clm},2));  % first take the mean of each trial, the mean accross trials  
                currErr(1, prev) = std(mean(allSesSuccessLick{cdIdx,clm},2))/sqrt(size(allSesSuccessLick{cdIdx,clm},1));  % standard error across alll trials
                
                numSamples = length(mean(allSesSuccessLick{cdIdx,clm},2));
                % collect data and parameters from current conditions for csv file
                csvOutput = [csvOutput; [ones(numSamples, 1), ...   % first column is subject ID
                                         repmat([preEV(prev), preOut(prot)], numSamples,1), ... % second and third columns are preEV and preOutcome                                 
                                         mean(allSesSuccessLick{cdIdx,clm},2)],...  % fourth column is data
                                                                              ];
                clear cdIdx numSamples
            end
            txtContent{prot} = ['\bf \color[rgb]{' num2str(colorCode(prot,:)) '} preOutcome = ' num2str(preOut(prot))];
            errorbar(1:3, currCrv, currErr, 'Color', colorCode(prot,:), 'LineWidth', 2)
            clear prev currCrv currErr 
        end
        % write the csv file and start from the second row
        csvwrite(csvFileName, [[1111 2222 3333 4444]; csvOutput], 0, 0);
        clear prot csvOutput csvFileName
        set(fg, 'XLim', [0 4], 'YLim', [0 1], 'XTick', 1:3, 'XTickLabel', {'3','6','9'})
        title([alignEvent{clm} ' [' num2str(eventTimeInterval{clm}(1)) ' ' num2str(eventTimeInterval{clm}(2)) ']'])
        xlabel('preEV')
        ylabel('mean lick signal')
        text(0.01, 0.45, txtContent)
        clear txtContent
    end
end
clear rw clm

% ....... Run ART manually ..............

% apply ANOVA on ART rank data and add the result to the subplots
for rw = 1:numRows
    for clm = 1:numColumns
        csvFileName = [csv_save_path monkeyName num2str(rw) num2str(clm) '.art.csv']; % full path plus file name for writing data of this frame in the csv file
        artData = csvread(csvFileName,1); % .art data -> [subjID preEV preOutcome Data, A(preEV) A(preOutcome) A(interaction), AR(preEV) AR(preOutcome) AR(interaction)] (column indices start from 0)
        % run parametric anova on A.R. data for preEV
        [~,tbl_preEV] = anovan(artData(:,8), artData(:,2:3),'model','interaction','varnames',{'preEV','preOutcome'},'display','off');
        txt = {['p-' tbl_preEV{2,1} '= ' num2str(tbl_preEV{2,7})]};
        % run parametric anova on A.R. data for preOutcome
        [~,tbl_preOut] = anovan(artData(:,9), artData(:,2:3),'model','interaction','varnames',{'preEV','preOutcome'},'display','off');
        txt = [txt; {['p-' tbl_preOut{3,1} '= ' num2str(tbl_preOut{3,7})]}];
        % run parametric anova on A.R. data for interaction
        [~,tbl_interact] = anovan(artData(:,10), artData(:,2:3),'model','interaction','varnames',{'preEV','preOutcome'},'display','off');
        txt = [txt; {['p-' tbl_interact{4,1} '= ' num2str(tbl_interact{4,7})]}];
        % add the final non-parametric anova test result to the current subplot
        subplot(numRows, numColumns, ((rw-1)*numColumns)+clm), text(1.2,0.2,txt);
        clear txt tbl_preEV tbl_preOut tbl_interact
    end
end

% add more info to the figure
subplot(numRows, numColumns,1)
text(-1,0.1,['\bf \fontsize{17} currEV = 3'],'Rotation',90)
subplot(numRows, numColumns,5)
text(-1,0.1,['\bf \fontsize{17} currEV = 6'],'Rotation',90)
subplot(numRows, numColumns,9)
text(-1,0.1,['\bf \fontsize{17} currEV = 9'],'Rotation',90)

%% prep data for GLM and run the glmfit command
allSesGLMpred = cell(length(allSesGroupedData), 1); % size = numSessions x 1; allSesGLMpred{ses,1} = [preOutcom preTotalReward preVar preEV currVar currEV]
allSesGLMdata = cell(length(allSesGroupedData), length(alignEvent)); % size = numSessions x numTimeIntervals; 

for ses = 1:length(allSesGroupedData) % for each session
    for gr = 1:length(allSesGroupedData{ses})  % for each group of trials
        allSesGLMpred{ses,1} = [ allSesGLMpred{ses,1}; [allSesGroupedData{ses}(gr).PreOutcome, allSesGroupedData{ses}(gr).PreTotalReward,  ...
                                                        allSesGroupedData{ses}(gr).PreVar, allSesGroupedData{ses}(gr).PreEV, ...
                                                        allSesGroupedData{ses}(gr).CurrVar, allSesGroupedData{ses}(gr).CurrEV] ];
       % for successful trials of this group, in each time interval, get the mean lick signal in that interval                                                    
        for tim = 1:length(alignEvent) % for each time interval
            eval([ 'allSesGLMdata{ses,tim} = [allSesGLMdata{ses,tim}; mean(allSesGroupedData{ses}(gr).SuccessLickInterval' num2str(tim) ',2)];' ]);           
        end
    end
end
clear ses gr tim

%% run GLM for each time interval
clmName = {[]};
for tim = 1:length(alignEvent) % for each time interval
    eval([ '[b' num2str(tim) ',dev' num2str(tim) ',stats' num2str(tim) '] = glmfit( cell2mat(allSesGLMpred), eps+cell2mat(allSesGLMdata(:,tim)), ''gamma'' );' ])
    clmName = {clmName{1,:}, ['Beta_time' num2str(tim)], ['pValue_time' num2str(tim)]};
end
clear tim
 clmName(1)=[];

for i = 1:length(clmName)/2
    eval([clmName{1+(i-1)*2} ' = b' num2str(i) ';'])
    eval([clmName{i*2} ' = stats' num2str(i) '.p;'])
end

RowNames = {'Constant';'Pre. outcome';'Pre. total reward'; 'Pre. Variance'; 'Pre. EV'; 'Curr. Variance'; 'Curr. EV'};
% display the result
T = table(Beta_time1, pValue_time1, Beta_time2, pValue_time2, Beta_time3, pValue_time3, Beta_time4, pValue_time4, 'RowNames', RowNames)
