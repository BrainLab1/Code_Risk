
% initiated on 12.05.2018

clear all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder bhv_data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
bhv_data_folder             = 'Z:\data\Risk\Behavior\';
         
Monkey = 'MacDuff';  % 'Moj' or 'Mac'
groupingCriteria       = 'expecter_reward & RewardVariance' ;  % $$ change the save_dir if you changed this
numConditions = 9; 
Fs = 1000;   % assume 1K sampling rate
num_ch = 96;

switch Monkey
    case 'MacDuff'
        eventTimeInterval = {[0.75 0.95]};
        alignEvent = {'cue'};  
        timeOfInt = [{[0.2 0.6]}, {[0.2 0.6]} , {[0.4 0.6]}];  % time in sec. based on Cue alignment
        freqOfInt = [{[10 20]}  , {[4 10]}    , {[4 10]}]; % frequency in hertz
    case 'Mojo'
        eventTimeInterval = {[0.4 0.8]};
        alignEvent = {'cue'};  
        timeOfInt = [{[0.2 0.4]}, {[0.4 0.6]}];  % time in sec. based on Cue alignment
        freqOfInt = [{[10 20]}  , {[4 10]}]; % frequency in hertz
end

if strcmp(Monkey(1:3),'Mac')
    clean_data_dir = 'Z:\data\Risk\Clean Data\';
    freq_data_dir = 'H:\Risk Task\MacDuff\';  % Bahareh PC
%     freq_data_dir = 'F:\Electrode-day power spectrum\Risk Task\MacDuff\';% Saeed PC
% % %     save_result_dir = 'X:\Report 15\GLM on electrode days\MacDuff\';
elseif strcmp(Monkey(1:3),'Moj')
    clean_data_dir = 'Z:\data\Risk\Clean Data\';
    freq_data_dir = 'H:\Risk Task\Mojo\';  % Bahareh PC
%     freq_data_dir = 'F:\Electrode-day power spectrum\Risk Task\Mojo\';% Saeed PC
% % % %     save_result_dir = 'X:\Report 15\GLM on electrode days\Mojo\';
end

ApplyLogTransformToPower = 1;
DiscriminabilityType = 'Cohen''s d';
save_dir = 'X:\Report 17\LFP vs Lick discriminability\Discriminability_[EV Var] grouping\';

%% Read out list of all the files related to this session
allFiles = dir(bhv_data_folder);    

% remove the '.' and '..' from the folders list and filter sessions for monkey names
idx = [];
for ses = 1:numel(allFiles) % for each session
    if (numel(allFiles(ses).name) < 3) || ~(strcmp(allFiles(ses).name(5:7), Monkey(1:3)))  
        idx = [idx; ses];
    end
end
allFiles(idx) = [];
clear ses idx
%%
allSesGroupedData = {};
for ses = 1:numel(allFiles) % for each session
    display(['Session ' num2str(ses)])

    % this will load new_cfg variable
    load ([bhv_data_folder allFiles(ses).name])
    
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
            % get the mean lick signal for each trial
            eval(['output(gr).SuccessLickInterval' num2str(tim) ' = mean(temp,2);'])
            eval(['output(gr).OrigCutSampleInterval' num2str(tim) ' = cutInterval;'])
            clear temp tr cutInterval
        end
        clear gr
        
        % ---- zscore mean lick signals in this session
        % get number of trials in each group
        eval(['sz = cellfun(@(x) length(x), {output.SuccessLickInterval' num2str(tim) '}'');'])
        % zscore mean lick across trials
        bb = struct2table(output); 
        eval(['bb.SuccessLickInterval' num2str(tim) ' = mat2cell(zscore(cell2mat({output.SuccessLickInterval' num2str(tim) '}'')), sz, 1);'])
        clear output
        output = table2struct(bb);
        clear bb sz
        
        for gr = 1:length(output)
            % collect lick signals of this group in this session's variable
            eval(['allSesSuccessLick(gr,tim) = {[allSesSuccessLick{gr,tim}; output(gr).SuccessLickInterval' num2str(tim) ']};'])
        end
        clear gr
    end
    allSesGroupedData{ses} = output;
    clear gr output eventTable new_cfg tmp_cfg temp tim
end
clear ses

%%  Get the discriminability index across conditions for the lich signal
lickDiscrMatrix = cell(size(allSesSuccessLick,2)); % length(lickDiscrMatrix) = number of time intervals
for tim = 1:length(lickDiscrMatrix)
   % initialize the discriminablity matrix
    lickDiscrMatrix{tim} = nan(size(allSesSuccessLick,1)); 
% % % %     lickDiscrMatrix{tim} = cell2mat( cellfun(@(x,y)   cell2mat(cellfun(@(x1,y)  (mean(x1)-mean(y))/sqrt((var(x1)/length(x1))+(var(y)/length(y)))  , repmat({x},size(allSesSuccessLick,1),1),  allSesSuccessLick(:,tim), 'UniformOutput', 0))'  , allSesSuccessLick(:,tim), 'UniformOutput', 0) );
    switch DiscriminabilityType 
        case 'Cohen''s d'
            lickDiscrMatrix{tim} = cell2mat( cellfun(@(x,y)   cell2mat(cellfun(@(x1,y)  (mean(x1)-mean(y))/sqrt( 0.5*((nanstd(x1)^2)+(nanstd(y)^2)) )  , repmat({x},size(allSesSuccessLick,1),1),  allSesSuccessLick(:,tim), 'UniformOutput', 0))'  , allSesSuccessLick(:,tim), 'UniformOutput', 0) );
    end
end
clear tim

%% read our relevant lfp file names
Array_pos = {'Anterior','Posterior'};% Array1 is anterior and Array2 is posterior
eval(['[Array1, Array2] = new_electrodepinout' Monkey '(''chan'');']) 
        
% read out all file names in the clean data folder
cleanDataFiles = dir(clean_data_dir);
cleanDataFiles(1:2) = [];
cleanDataFiles = {cleanDataFiles.name}';

% select file names related to the selected monkey
monkeyCleanData = cleanDataFiles(~cellfun(@isempty, strfind(cleanDataFiles, Monkey)));   

freqDataFiles = dir(freq_data_dir);
freqDataFiles(1:2) = [];
freqDataFiles = {freqDataFiles.name}';

% find indices for frequency and time of interest
% read out one of the file names related to 96 channels of this session
tmpIdx1 = strfind(monkeyCleanData{1},'M');
tmpIdx2 = strfind(monkeyCleanData{1},'(')-1;
sesChFileNames = freqDataFiles( ~cellfun( @isempty, strfind(freqDataFiles, monkeyCleanData{1}(tmpIdx1:tmpIdx2)) ) );
clear tmpIdx1 tmpIdx2
load([freq_data_dir sesChFileNames{1}])
if ~exist('allFreqTimIdx','var')
    % find the indices for desired freq 
    [~,a] = cellfun(@(x) min(abs(x(1)-freq.freq)) , freqOfInt);
    fIdx(:,1) = a';
    clear a
    [~,a] = cellfun(@(x) min(abs(x(2)-freq.freq)) , freqOfInt);
    fIdx(:,2) = a';
    clear a
    % find the indices for desired time 
    [~,a] = cellfun(@(x) min(abs(x(1)-freq.time)) , timeOfInt);
    tIdx(:,1) = a';
    clear a
    [~,a] = cellfun(@(x) min(abs(x(2)-freq.time)) , timeOfInt);
    tIdx(:,2) = a';
    clear a
    % get all combinations of freq and time of interest intervals -> allFreTimInterest(i,:) = [freInteresStart freInteresEnd timInterestStart timInterestEnd]
    allFreqTimIdx = [fIdx,tIdx];
            
    allROIinfo = {};
    for tf=1:size(allFreqTimIdx,1)
        allROIinfo{tf} = ['freq[' num2str(round(freq.freq(allFreqTimIdx(tf,1)))) ',' num2str(round(freq.freq(allFreqTimIdx(tf,2)))), ...
                          '] Cue[' num2str(freq.time(allFreqTimIdx(tf,3))) ',' num2str(freq.time(allFreqTimIdx(tf,4))) ']' ];
    end
    clear tIdx fIdx tf
end
clear freq

% initialize the variable!!
allSesMeanPower1 = cell(numConditions, size(allFreqTimIdx,1)); % size= [numGroups (numFreInterest x numTimeInterest)]
allSesMeanPower2 = cell(numConditions, size(allFreqTimIdx,1)); % size= [numGroups (numFreInterest x numTimeInterest)]
arrayIdx = nan(length(monkeyCleanData),num_ch);

for ses = 1:length(monkeyCleanData) 
    % read out file names related to 96 channels of this session
    tmpIdx1 = strfind(monkeyCleanData{ses},'M');
    tmpIdx2 = strfind(monkeyCleanData{ses},'(')-1;
    sesChFileNames = freqDataFiles( ~cellfun( @isempty, strfind(freqDataFiles, monkeyCleanData{ses}(tmpIdx1:tmpIdx2)) ) );
    clear tmpIdx1 tmpIdx2
    
    % sort the file names based on the channel lable index (i.e., X in the 'chX' of the freq. file name !!)
    a1 = cellfun(@(x)    str2double(x(max(strfind(x,'ch'))+2 : min(strfind(x,'.'))-1)), sesChFileNames);
    [~,sortedIdx]=sort(a1);
    sesChFileNames = sesChFileNames(sortedIdx);
    clear a1 sortedIdx
    
    for ch = 1:length(sesChFileNames)
        display(['Ses ' num2str(ses) ' , Ch ' num2str(ch)])
        
        % load the electrode power spectrum; this will load the freq variable ... 
        load([freq_data_dir sesChFileNames{ch}])
        % find which array this channel belongs to
        if ismember(ch, Array1)
            arrayIdx(ses,ch) = 1;
        else if ismember(ch, Array2)
                arrayIdx(ses,ch) = 2;
            else
                display('There is an error! This channels belong to neither of the arrays!!!')
                break
            end
        end
        
        % change order of freq.powspctrm dimensions and change it to [numFreqBands, numTimeStamps, numTrials]
        M0 = permute(squeeze(freq.powspctrm(:,1,:,:)), [2,3,1]);
        % turn the M0 matrix to a 2D matrix of size [numFreqBands (numTrials*numTimaStamps)]
        M = reshape(M0, size(M0,1), []);
        % zscore each frequency band separately and exclude NaN powers (i.e., trials in which electrode was outlier)
        Z_M = (M - repmat(nanmean(M,2),1,size(M,2)))./repmat(nanstd(M,1,2),1,size(M,2));
        % reshape the zscored data into a 3D matrix of size [numFreqBands, numTimeStamps, numTrials]
        Z_freq_powspctrm = reshape(Z_M,size(M0));
        clear M0 M Z_M

        % apply log transformation if asked
       switch ApplyLogTransformToPower
           case 1
               Z_freq_powspctrm = log(Z_freq_powspctrm);
       end
           
        % read out the event structure
        event = recursive_read(freq, 'event');  % this is the same event structure which were saved in the Clean Data folder for this session
        
        % group trials based on the specified criterion
        [output] = GroupTrials(event, groupingCriteria);

        for freqtim = 1:size(allFreqTimIdx)
           aa = cellfun(@(x) nanmean(reshape(x,[],size(x,3)),1)'  , ...
                   cellfun(@(x) Z_freq_powspctrm( allFreqTimIdx(freqtim,1):allFreqTimIdx(freqtim,2), allFreqTimIdx(freqtim,3):allFreqTimIdx(freqtim,4), x), {output.TrialIdx}', 'UniformOutput', 0) , ...
                   'UniformOutput', 0) ;
           switch arrayIdx(ses,ch)
               case 1
                   allSesMeanPower1(:,freqtim) = cellfun(@(x,y) [y;x] , allSesMeanPower1(:,freqtim), aa, 'UniformOutput', 0);
               case 2
                   allSesMeanPower2(:,freqtim) = cellfun(@(x,y) [y;x] , allSesMeanPower2(:,freqtim), aa, 'UniformOutput', 0);
           end
           clear aa
        end
        
        clear freqtim output event Z_freq_powspctrm freq
    
    end  % end of ch=1:...
    
    clear sesChFileNames
end

%%  Get the discriminability index across conditions for LFP power
lfpDiscrMatrix1 = cell(size(allFreqTimIdx,1)); % length(lfpkDiscrMatrix1) = number of time intervals
lfpDiscrMatrix2 = cell(size(allFreqTimIdx,1)); % length(lfpkDiscrMatrix2) = number of time intervals
for freqtim = 1:length(lfpDiscrMatrix1)
   % initialize the discriminablity matrix
    lfpDiscrMatrix1{freqtim} = nan(size(allSesMeanPower1,1)); 
    lfpDiscrMatrix2{freqtim} = nan(size(allSesMeanPower2,1)); 

% % % %     lfpDiscrMatrix1{freqtim} = cell2mat( cellfun(@(x,y)   cell2mat(cellfun(@(x1,y)  (nanmean(x1)-nanmean(y))/sqrt((nanvar(x1)/sum(~isnan(x1)))+(nanvar(y)/sum(~isnan(y))))  , repmat({x},size(allSesMeanPower1,1),1),  allSesMeanPower1(:,freqtim), 'UniformOutput', 0))'  , allSesMeanPower1(:,freqtim), 'UniformOutput', 0) );
% % % %     lfpDiscrMatrix2{freqtim} = cell2mat( cellfun(@(x,y)   cell2mat(cellfun(@(x1,y)  (nanmean(x1)-nanmean(y))/sqrt((nanvar(x1)/sum(~isnan(x1)))+(nanvar(y)/sum(~isnan(y))))  , repmat({x},size(allSesMeanPower2,1),1),  allSesMeanPower2(:,freqtim), 'UniformOutput', 0))'  , allSesMeanPower2(:,freqtim), 'UniformOutput', 0) );
    switch DiscriminabilityType 
        case 'Cohen''s d'
            lfpDiscrMatrix1{freqtim} = cell2mat( cellfun(@(x,y)   cell2mat(cellfun(@(x1,y)  (nanmean(x1)-nanmean(y))/sqrt( 0.5*((nanstd(x1)^2)+(nanstd(y)^2)) )  , repmat({x},size(allSesMeanPower1,1),1),  allSesMeanPower1(:,freqtim), 'UniformOutput', 0))'  , allSesMeanPower1(:,freqtim), 'UniformOutput', 0) );
            lfpDiscrMatrix2{freqtim} = cell2mat( cellfun(@(x,y)   cell2mat(cellfun(@(x1,y)  (nanmean(x1)-nanmean(y))/sqrt( 0.5*((nanstd(x1)^2)+(nanstd(y)^2)) )  , repmat({x},size(allSesMeanPower2,1),1),  allSesMeanPower2(:,freqtim), 'UniformOutput', 0))'  , allSesMeanPower2(:,freqtim), 'UniformOutput', 0) );
    end
end
clear tim

%% plot lick vs lfp discriminability

fg1 = figure('Name', [Monkey ', dlPFC'])
numRows = 2;
numColumns = 2;
plotCounter = 0;
for lickTimIntrvl = 1:length(lickDiscrMatrix)
    for lfpFreqTim = 1:length(lfpDiscrMatrix1)
       plotCounter = plotCounter+1;
       subplot(numRows, numColumns, plotCounter), hold on, box on
       titl = [ {['lick ' alignEvent{lickTimIntrvl} '[' num2str(eventTimeInterval{lickTimIntrvl}(1)) ' , ' num2str(eventTimeInterval{lickTimIntrvl}(2)) ']' ]} , {['LFP: ' allROIinfo{lfpFreqTim}]} ];
       switch ApplyLogTransformToPower
           case 1
               titl = [{'power log transformed'} , titl];
       end
       title(titl)
       xlabel('lick discriminability index')
       ylabel('lfp discriminability index')
       
       for cd1 = 1 : size(lickDiscrMatrix{lickTimIntrvl},1)-1
           for cd2 = cd1+1 : size(lickDiscrMatrix{lickTimIntrvl},1)
               scatter(abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) , abs(lfpDiscrMatrix1{lfpFreqTim}(cd1,cd2)), 'Marker', 'o', 'MarkerEdgeColor', 'k')
               scatter(abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) , abs(lfpDiscrMatrix2{lfpFreqTim}(cd1,cd2)), 'Marker', 'o', 'MarkerEdgeColor', 'r')
           end
           clear cd2
       end
       clear cd1 titl
       plot([min([xlim ylim]) max([xlim ylim])],[min([xlim ylim]) max([xlim ylim])],'--', 'Color', 'b')
       axis square 
        
    end
end
close(fg1)

%%
colorCode = [0 1 1
             0 0 1
             0 0 0
             ...
             1 0.6 0.6    
             0.99 0.2 0.2
             0.6 0 0
             ...
             0.6 1 0.6
             0 1 0
             0.1 0.5 0.1];

% get the parameter values for each condition
grpParaVal = {allSesGroupedData{1,1}.Value}';
% grpParaValString = cellfun(@(x) ['[' num2str(x(1)) ',' num2str(x(2)) ']'], grpParaVal, 'UniformOutput', 0);

fg2 = figure('Name', [Monkey ', 2 areas'])
fg2_dlpfc = figure('Name', [Monkey ', dlPFC'], 'Position', [2341 117 928 845])
fg2_7a = figure('Name', [Monkey ', 7a'], 'Position', [2341 117 928 845])

numRows = 2;
numColumns = 2;
plotCounter = 0;
for lickTimIntrvl = 1:length(lickDiscrMatrix)
    for lfpFreqTim = 1:length(lfpDiscrMatrix1)
       plotCounter = plotCounter+1;
       % scatter plot 2 areas
       figure(fg2)
       subplot(numRows, numColumns, plotCounter), hold on, box on
       titl = [ {['lick ' alignEvent{lickTimIntrvl} '[' num2str(eventTimeInterval{lickTimIntrvl}(1)) ' , ' num2str(eventTimeInterval{lickTimIntrvl}(2)) ']' ]} , {['LFP: ' allROIinfo{lfpFreqTim}]} ];
       switch ApplyLogTransformToPower
           case 1
               titl = [{'power log transformed'} , titl];
       end
       title(titl)
       xlabel('lick discriminability index')
       ylabel('lfp discriminability index')
       xTickLabel = [];
       pairedCdDiscrimIdxLick = [];
       pairedCdDiscrimIdxLFPArr1 = [];
       pairedCdDiscrimIdxLFPArr2 = [];
       pairedCdIdx = []; 
       for cd1 = 1 : size(lickDiscrMatrix{lickTimIntrvl},1)-1
           for cd2 = cd1+1 : size(lickDiscrMatrix{lickTimIntrvl},1)
%                xTickLabel = [ xTickLabel, {[grpParaValString{cd1} ' ' grpParaValString{cd2}]} ];
               pairedCdDiscrimIdxLick = [pairedCdDiscrimIdxLick, abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2))];
               pairedCdDiscrimIdxLFPArr1  = [pairedCdDiscrimIdxLFPArr1 , abs(lfpDiscrMatrix1{lfpFreqTim}(cd1,cd2))];
               pairedCdDiscrimIdxLFPArr2  = [pairedCdDiscrimIdxLFPArr2 , abs(lfpDiscrMatrix2{lfpFreqTim}(cd1,cd2))];
               pairedCdIdx = [pairedCdIdx; [cd1 cd2]];
              scatter(abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) , abs(lfpDiscrMatrix1{lfpFreqTim}(cd1,cd2)), 'Marker', 'o', 'MarkerEdgeColor', 'k')
              scatter(abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) , abs(lfpDiscrMatrix2{lfpFreqTim}(cd1,cd2)), 'Marker', 'o', 'MarkerEdgeColor', 'r')
              legend('dlPFC','7a')
           end
           clear cd2
       end
       clear cd1 titl        
       
       % scatter plot dlPFC
       figure(fg2_dlpfc)
       subplot(numRows, numColumns, plotCounter), hold on, box on
       titl = [ {['lick ' alignEvent{lickTimIntrvl} '[' num2str(eventTimeInterval{lickTimIntrvl}(1)) ' , ' num2str(eventTimeInterval{lickTimIntrvl}(2)) ']' ]} , {['LFP: ' allROIinfo{lfpFreqTim}]} ];
       switch ApplyLogTransformToPower
           case 1
               titl = [{'power log transformed'} , titl];
       end
       title(titl)
       xlabel('lick discriminability index')
       ylabel('lfp discriminability index')
       xVal = [];
       yVal = [];
       for cd1 = 1 : size(lickDiscrMatrix{lickTimIntrvl},1)-1
           for cd2 = cd1+1 : size(lickDiscrMatrix{lickTimIntrvl},1)
               % plot dlPFC
               plot(abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) , abs(lfpDiscrMatrix1{lfpFreqTim}(cd1,cd2)), 'Marker', 'o', ...
                   'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerEdgeColor', colorCode(cd1,:), 'MarkerFaceColor', colorCode(cd2,:))
               xVal = [xVal; abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) ];
               yVal = [yVal; abs(lfpDiscrMatrix1{lfpFreqTim}(cd1,cd2))];
           end
           clear cd2
       end
       [b,bint] = regress(yVal,[ones(size(xVal)) xVal]);
       plot(xlim, b(1)+(b(2)*xlim), 'k');
       text(mean(xlim), mean(ylim), [{['Slope ' num2str(b(2))]},{['conf. interval [' num2str(bint(2,1)) ',' num2str(bint(2,2)) ']']}])
       clear cd1 titl xVal yVal b bint    
       
       % scatter plot 7a
       figure(fg2_7a)
       subplot(numRows, numColumns, plotCounter), hold on, box on
       titl = [ {['lick ' alignEvent{lickTimIntrvl} '[' num2str(eventTimeInterval{lickTimIntrvl}(1)) ' , ' num2str(eventTimeInterval{lickTimIntrvl}(2)) ']' ]} , {['LFP: ' allROIinfo{lfpFreqTim}]} ];
       switch ApplyLogTransformToPower
           case 1
               titl = [{'power log transformed'} , titl];
       end
       title(titl)
       xlabel('lick discriminability index')
       ylabel('lfp discriminability index')
       xVal = [];
       yVal = [];
       for cd1 = 1 : size(lickDiscrMatrix{lickTimIntrvl},1)-1
           for cd2 = cd1+1 : size(lickDiscrMatrix{lickTimIntrvl},1)
               % plot 7a
               plot(abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) , abs(lfpDiscrMatrix2{lfpFreqTim}(cd1,cd2)), 'Marker', 'o', ...
                   'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerEdgeColor', colorCode(cd1,:), 'MarkerFaceColor', colorCode(cd2,:))
               xVal = [xVal; abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) ];
               yVal = [yVal; abs(lfpDiscrMatrix2{lfpFreqTim}(cd1,cd2))];
           end
           clear cd2
       end
       [b,bint] = regress(yVal,[ones(size(xVal)) xVal]);
       plot(xlim, b(1)+(b(2)*xlim), 'k');
       text(mean(xlim), mean(ylim), [{['Slope ' num2str(b(2))]},{['conf. interval [' num2str(bint(2,1)) ',' num2str(bint(2,2)) ']']}])
       clear cd1 titl b bint      
       
       
       % line plot
       if 0
        plot(1:size(pairedCdIdx,1) , [sign(diff(pairedCdDiscrimIdxLick)) 0]  , 'Color', 'r')
        plot(1:size(pairedCdIdx,1) , [sign(diff(pairedCdDiscrimIdxLFPArr1)) 0], 'Color', 'c')
        plot(1:size(pairedCdIdx,1) , [sign(diff(pairedCdDiscrimIdxLFPArr2)) 0], 'Color', 'b')
       
        set(gca, 'XTick', 1:size(pairedCdIdx,1), 'XTickLabel', xTickLabel)
        xtickangle(gca,90)
        legend('Lick','dlPFC','7a')
       end
       
    end
end
close(fg2)

saveas(fg2_dlpfc, [save_dir Monkey '_Scatter_LFPDiscVsLick_dlPFC'], 'fig')
saveas(fg2_dlpfc, [save_dir Monkey '_Scatter_LFPDiscVsLick_dlPFC'], 'bmp')

saveas(fg2_7a, [save_dir Monkey '_Scatter_LFPDiscVsLick_7a'], 'fig')
saveas(fg2_7a, [save_dir Monkey '_Scatter_LFPDiscVsLick_7a'], 'bmp')

close(fg2_dlpfc), clear fg2_dlpfc
close(fg2_7a), clear fg2_7a
