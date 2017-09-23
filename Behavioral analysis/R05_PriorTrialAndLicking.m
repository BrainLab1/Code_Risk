
% this script plots and analyzes effect of the prior trial EV and Outcome on
% the licking probability of the current trial, in four time intervals 
% groups trials based on 'PreTrlVar & PreTrlOutcome & CurrTrlEV'
% runs GLM with following factors: {'Constant';'Pre. outcome';'Pre. Variance'; 'Pre. EV'; 'Curr. Variance'; 'Curr. EV'}
% GLM regressors are not zscored!
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
monkeyName = 'Mac';  % 'Moj' or 'Mac'
csv_save_path           = ['D:\Report05_ART for licking data\' monkeyName '\'];
groupingCriteria       = 'PreTrlVar & PreTrlOutcome & CurrTrlEV' ; 
Fs = 1000;   % assume 1K sampling rate

% each row of the 'allConditions' contains following values -> 
%                [PreTrlVar   PreTrlOutcome  CurrTrlEV] 
allConditions = [0  0   3   %  -
                 0  0   6   % | 
                 0  0   9   %  -
                 ...
                 1  -1  3   %  -
                 1  -1  6   % | 
                 1  -1  9   %  -
                 1  1   3   %  -
                 1  1   6   % | 
                 1  1   9   %  -
                 ...
                 4  -1  3   %  -
                 4  -1  6   % | 
                 4  -1  9   %  -
                 4  1   3   %  -
                 4  1   6   % | 
                 4  1   9]; %  -
                     
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
            output(gr).PreEV = eventTable.expected_reward(output(gr).SuccessTrlIdx - 1);
            % get the previous trial outcome 
            output(gr).PreOutcome = repmat(output(gr).Value(2), length(output(gr).TrialIdx),1);
            % get the previous trial expected value
            output(gr).PreVar = repmat(output(gr).Value(1), length(output(gr).TrialIdx),1);
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
preVar = [1;4];
preOut = [-1;1];  % previouse-trila outcome [lose neutral win]
pltCounter = 0;

for rw = 1:numRows
    for clm = 1:numColumns
        pltCounter = pltCounter+1;
        fg = subplot(numRows, numColumns, pltCounter); hold on, box on
        txtContent = {};  % this is the text which will show the legend
        csvOutput = [];  % this is the matrix which will be written in the csv file
        csvFileName = [csv_save_path monkeyName num2str(rw) num2str(clm) '.csv']; % full path plus file name for writing data of this frame in the csv file
        
        % get condition [0 0 plotRowCurrEV(rw)]
        cdIdx0 = find( (allConditions(:,3) == plotRowCurrEV(rw)) .* (allConditions(:,2) == 0) .* (allConditions(:,1) == 0) );
        currCrv0 = mean(mean(allSesSuccessLick{cdIdx0,clm},2));  % first take the mean of each trial, the mean accross trials  
        currErr0 = std(mean(allSesSuccessLick{cdIdx0,clm},2))/sqrt(size(allSesSuccessLick{cdIdx0,clm},1));  % standard error across alll trials
        numSamples0 = length(mean(allSesSuccessLick{cdIdx0,clm},2));

        errorbar(1, currCrv0, currErr0, '.', 'Color', colorCode(3,:), 'LineWidth', 2)
        txtContent{3} = ['\bf \color[rgb]{' num2str(colorCode(3,:)) '} preOutcome = 0' ];
        clear cdIdx0 numSamples0 currCrv0 currErr0                                                                
        
        for prot = 1:length(preOut)  %for each preOut
            currCrv = [];  % current curve!
            currErr = [];  % current error
            % form the curve with thrwee points 
            for prvar = 1:length(preVar) % for each preEV
                % find conditions to be ploted in this row
                cdIdx = find( (allConditions(:,3) == plotRowCurrEV(rw)) .* (allConditions(:,2) == preOut(prot)) .* (allConditions(:,1) == preVar(prvar)) );
                currCrv(1, prvar) = mean(mean(allSesSuccessLick{cdIdx,clm},2));  % first take the mean of each trial, the mean accross trials  
                currErr(1, prvar) = std(mean(allSesSuccessLick{cdIdx,clm},2))/sqrt(size(allSesSuccessLick{cdIdx,clm},1));  % standard error across alll trials
                
                numSamples = length(mean(allSesSuccessLick{cdIdx,clm},2));
                % collect data and parameters from current conditions for csv file
                csvOutput = [csvOutput; [ones(numSamples, 1), ...   % first column is subject ID
                                         repmat([preVar(prvar), preOut(prot)], numSamples,1), ... % second and third columns are preVar and preOutcome                                 
                                         mean(allSesSuccessLick{cdIdx,clm},2)],...  % fourth column is data
                                                                              ];
                clear cdIdx numSamples
            end
            txtContent{prot} = ['\bf \color[rgb]{' num2str(colorCode(prot,:)) '} preOutcome = ' num2str(preOut(prot))];
            errorbar(2:length(currCrv)+1, currCrv, currErr, 'Color', colorCode(prot,:), 'LineWidth', 2)
            clear prvar currCrv currErr 
        end
        % write the csv file and start from the second row
        csvwrite(csvFileName, [[1111 2222 3333 4444]; csvOutput], 0, 0);
        clear prot csvOutput csvFileName
        set(fg, 'XLim', [0 4], 'YLim', [0 1], 'XTick', 1:3, 'XTickLabel', {'0','1','4'})
        title([alignEvent{clm} ' [' num2str(eventTimeInterval{clm}(1)) ' ' num2str(eventTimeInterval{clm}(2)) ']'])
        xlabel('preVar')
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
        txt = {'non-parametric ANOVA';'-----------------'};
        csvFileName = [csv_save_path monkeyName num2str(rw) num2str(clm) '.art.csv']; % full path plus file name for writing data of this frame in the csv file
        artData = csvread(csvFileName,1); % .art data -> [subjID preVar preOutcome Data, A(preVar) A(preOutcome) A(interaction), AR(preVar) AR(preOutcome) AR(interaction)] (column indices start from 0)
        % run parametric anova on A.R. data for preVar
        [~,tbl_preVar] = anovan(artData(:,8), artData(:,2:3),'model','interaction','varnames',{'preVar','preOutcome'},'display','off');
        txt = [txt; {['p-' tbl_preVar{2,1} '= ' num2str(tbl_preVar{2,7})]}];
        % run parametric anova on A.R. data for preOutcome
        [~,tbl_preOut] = anovan(artData(:,9), artData(:,2:3),'model','interaction','varnames',{'preVar','preOutcome'},'display','off');
        txt = [txt; {['p-' tbl_preOut{3,1} '= ' num2str(tbl_preOut{3,7})]}];
        % run parametric anova on A.R. data for interaction
        [~,tbl_interact] = anovan(artData(:,10), artData(:,2:3),'model','interaction','varnames',{'preVar','preOutcome'},'display','off');
        txt = [txt; {['p-' tbl_interact{4,1} '= ' num2str(tbl_interact{4,7})]}];
        % add the final non-parametric anova test result to the current subplot
        subplot(numRows, numColumns, ((rw-1)*numColumns)+clm), text(1.5,0.25,txt);
        clear txt tbl_preVar tbl_preOut tbl_interact
    end
end

% add more info to the figure
subplot(numRows, numColumns,1)
text(-1,0.1,['\bf \fontsize{17} currEV = 3'],'Rotation',90)
subplot(numRows, numColumns,5)
text(-1,0.1,['\bf \fontsize{17} currEV = 6'],'Rotation',90)
subplot(numRows, numColumns,9)
text(-1,0.1,['\bf \fontsize{17} currEV = 9'],'Rotation',90)
subplot(numRows, numColumns,8); 
text(4.5,2.5,['\bf \fontsize{17} non-para. ANOVA performed only on preVar 1 & 4, preOutcome win and lose'],'Rotation',-90)

%% prep data for GLM and run the glmfit command
allSesGLMpred = cell(length(allSesGroupedData), 1); % size = numSessions x 1; allSesGLMpred{ses,1} = [preOutcom preVar preEV currVar currEV]
allSesGLMdata = cell(length(allSesGroupedData), length(alignEvent)); % size = numSessions x numTimeIntervals; 

allSesGLMpred2 = cell(length(allSesGroupedData), 1); % size = numSessions x 1; allSesGLMpred{ses,1} = [preTotalReward currVar currEV]

for ses = 1:length(allSesGroupedData) % for each session
    for gr = 1:length(allSesGroupedData{ses})  % for each group of trials
        allSesGLMpred{ses,1} = [ allSesGLMpred{ses,1}; [allSesGroupedData{ses}(gr).PreOutcome, ...
                                                        allSesGroupedData{ses}(gr).PreVar, allSesGroupedData{ses}(gr).PreEV, ...
                                                        allSesGroupedData{ses}(gr).CurrVar, allSesGroupedData{ses}(gr).CurrEV] ];
        allSesGLMpred2{ses,1} = [ allSesGLMpred2{ses,1}; [allSesGroupedData{ses}(gr).PreTotalReward, ...
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
%     eval([ '[b' num2str(tim) ',dev' num2str(tim) ',stats' num2str(tim) '] = glmfit( cell2mat(allSesGLMpred), eps+cell2mat(allSesGLMdata(:,tim)), ''gamma'' );' ])
    eval([ '[b' num2str(tim) ',dev' num2str(tim) ',stats' num2str(tim) '] = glmfit( cell2mat(allSesGLMpred), cell2mat(allSesGLMdata(:,tim)), ''normal'' );' ])
% % % %     eval([ '[bb' num2str(tim) ',devv' num2str(tim) ',statss' num2str(tim) '] = glmfit( cell2mat(allSesGLMpred2), eps+cell2mat(allSesGLMdata(:,tim)), ''gamma'' );' ])
    clmName = {clmName{1,:}, ['Beta_time' num2str(tim)], ['pValue_time' num2str(tim)], ['pUnderThrsh' num2str(tim)]};
end
clear tim
clmName(1)=[];

for i = 1:length(clmName)/3
    eval([clmName{1+(i-1)*3} ' = b' num2str(i) ';'])
    eval([clmName{2+(i-1)*3} ' = stats' num2str(i) '.p;'])
    % find the false discovery rate threshold
    eval(['[pthr,pcor] = fdr(pValue_time' num2str(i) ');'])
    eval([clmName{i*3} ' = pValue_time' num2str(i) ' <= pthr;']);
    
% % % %     eval(['X' clmName{1+(i-1)*2} ' = bb' num2str(i) ';'])
% % % %     eval(['X' clmName{i*2} ' = statss' num2str(i) '.p;'])
end

RowNames = {'Constant';'Pre. outcome';'Pre. Variance'; 'Pre. EV'; 'Curr. Variance'; 'Curr. EV'};

% display the result
T = table(Beta_time1, pValue_time1, pUnderThrsh1, Beta_time2, pValue_time2, pUnderThrsh2,...
          Beta_time3, pValue_time3, pUnderThrsh3, Beta_time4, pValue_time4, pUnderThrsh4, 'RowNames', RowNames)


% write the GLM result in the xls file
% % % % xlswrite(['D:\Dropbox\Shared_Jackie Project\Report 5\' monkeyName 'GLM'], T{1:end,1:end}, 'B2:I7')


% % % % RowNames2 = {'Constant';'Pre. TotalReward';'Curr. Variance'; 'Curr. EV'};
% % % % % display the result
% % % % T2 = table(XBeta_time1, XpValue_time1, XBeta_time2, XpValue_time2, XBeta_time3, XpValue_time3, XBeta_time4, XpValue_time4, 'RowNames', RowNames2)
