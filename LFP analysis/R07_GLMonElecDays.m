

% This script loads LFP and event data from the Clean Data folder; groups trials based on 'expecter_reward & RewardVariance';  
% For each electrode in each condition, determines in which trials the electrode had bad signal and has to be removed
% from the analysis; it removes all electrodes from the analysis for trials with wierd (outlier!) event timing;
% In each session, we applied GLM on those channels which had enough repetition-per-group for all groups 
% GLM input data was mean power in the time and freq. interval of interest, measured separately in each trial
% We zscored the parameter values which were used for grouping trials and added a new field .ZValue to the output. These
% values will serve as GLM predictors
 


clear
close all
clc

%% set some parameters
Monkey            = 'Mojo';
groupingCriteria  = 'expecter_reward & RewardVariance';
num_ch = 96;
minElecRepPerCnd = 10;
timeOfInt = [{[0 0.4]}];  % time in sec. based on Cue alignment
freqOfInt = [{[15 20]},{[25 30]}]; % frequency in hertz

%%
if strcmp(Monkey,'MacDuff')
    clean_data_dir = 'Z:\data\Risk\Clean Data\';
%     freq_data_dir = 'H:\MacDuff\'; Bahareh PC
    freq_data_dir = 'F:\Electrode-day power spectrum\MacDuff\';% Saeed PC
    save_result_dir = 'X:\Report 7\Power analysis\MacDuff\';
elseif strcmp(Monkey,'Mojo')
    clean_data_dir = 'Z:\data\Risk\Clean Data\';
%     freq_data_dir = 'H:\Mojo\';% Bahareh PC
    freq_data_dir = 'F:\Electrode-day power spectrum\Mojo\';% Saeed PC
    save_result_dir = 'X:\Report 7\Power analysis\Mojo\';
end

%% read our relevant file names
Array_pos = {'Anterior','Posterior'};% Array1 is anterior and Array2 is posterior

% read out all file names in the clean data folder
cleanDataFiles = dir(clean_data_dir);
cleanDataFiles(1:2) = [];
cleanDataFiles = {cleanDataFiles.name}';

% select file names related to the selected monkey
monkeyCleanData = cleanDataFiles(~cellfun(@isempty, strfind(cleanDataFiles, Monkey)));   

%% find bad trial numbers and number of good trials for each electrode-group
for ses = 1:length(monkeyCleanData) 
    ses
    % load the current session data; this will load a variable named data which include both LFP signals and event info 
    load([clean_data_dir monkeyCleanData{ses}])
    
    % read out the event structure
    event = recursive_read(data,'event');
    
    % group trials based on the specified criterion
    [output] = GroupTrials(event, groupingCriteria);
    
    % convet event structure to table
    event = struct2table(event);
    
    % initialize the variable if doesn't exist !! this variable tells us in each condition, for which trials the electrode  
    % must be removed; it contains trial indices
    if ~exist('elecCndBadTrials','var') 
        elecCndBadTrials = cell(size(output,1), 96, length(monkeyCleanData)); % size: [numGroups x numElectrodes x numSessions]
    end
    
    % this variable shows number of valid trials per condition  for each electrode
    if ~exist('elecCndNumGoodTrials','var') 
        elecCndNumGoodTrials = nan(size(output,1), 96, length(monkeyCleanData) ); % size: [numGroups x numElectrodes x numSessions]
    end
    
    % fill in the elecCndBadTrials cells with trial numbers
    for gr=1:size(output,1)
        for str = 1:length(output(gr).TrialIdx) % for each trials within this group
            trIdx = output(gr).TrialIdx(str);
            
            % pool bad_lfp_channels and bad_amplitude_channels
            allBadCh = unique( [event.bad_lfp_channels{trIdx,1}'; event.bad_amplitude_channels{trIdx,1}] );
            
            % if the trial has to be removed because of event timing outliers, add all 96 channels to the allBadCh!
            if (event.OutCueOnOffZ2(trIdx) + event.OutAqToRwdZ2(trIdx))
                allBadCh = unique([allBadCh; [1:96]']);
            end
            
            % add this trial index to for the bad channels in this trial!!!
            elecCndBadTrials(gr,allBadCh,ses) =  cellfun( @(x) {[elecCndBadTrials{gr,x,ses}; trIdx]}, mat2cell(allBadCh,repmat(1,length(allBadCh),1),1) );
            clear trIdx allBadCh
        end
        % find number of good trials in this condition for each electrode 
        elecCndNumGoodTrials(gr,:,ses) = length(output(gr).TrialIdx) - cellfun(@length, elecCndBadTrials(gr,:,ses));
    end
    clear gr data event output
   
end

%% 
freqDataFiles = dir(freq_data_dir);
freqDataFiles(1:2) = [];
freqDataFiles = {freqDataFiles.name}';

% initialize the variable!!
allSesMeanPower = cell(length(monkeyCleanData), 9, length(freqOfInt), length(timeOfInt)); % size= [numSesssions x numGroups x numFreInterest x numTimeInterest]
glmB = cell(length(monkeyCleanData) ,96, 2);
glmStats = cell(length(monkeyCleanData) ,96, 2);

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
        
        % find which array this channel belongs to
        eval(['[Array1, Array2] = electrodepinout' Monkey '(freq.label{1}(1:4));']) 
        switch freq.label{1}(1:4)
            case 'chan'
                if ismember(ch, Array1)
                    arrayIdx(ch) = 1;
                else if ismember(ch, Array2)
                        arrayIdx(ch) = 2;
                    else
                        display('There is an error! This channels belong to neither of the arrays!!!')
                        break
                    end
                end
                
            case 'elec'
                k = strfind(freq.label{1},'-');
                if ismember(str2num(freq.label{1}(k+1:end)), Array1)
                    arrayIdx(ch) = 1;
                else if ismember(str2num(freq.label{1}(k+1:end)), Array2)
                        arrayIdx(ch) = 2;
                    else
                        udisplay('There is an error! This channels belong to neither of the arrays!!!')
                        break
                    end
                end
                clear k
        end
        
        % if electrode didn't have enough repetition-per-condition in all conditions, exclude the electrode from analysis
        if sum(elecCndNumGoodTrials(:,ch,ses) >= minElecRepPerCnd) ~= size(elecCndNumGoodTrials,1)
            continue
        end
        
        % load the electrode power spectrum if it has passed the above criterion; this will load the freq variable ... 
        load([freq_data_dir sesChFileNames{ch}])
        
        % find indices for frequency and time of interest
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
            allFreqTimIdx = combvec(fIdx',tIdx')';
            clear tIdx fIdx
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
        
        % read out the event structure
        event = recursive_read(freq, 'event');  % this is the same event structure which were saved in the Clean Data folder for this session
        
        % group trials based on the specified criterion
        [output] = GroupTrials(event, groupingCriteria);
        
        % zscore GLM predictors, get them from output and add a new field .ZValue to the output
        aa = cell2mat({output.Value}');
        [C1,~,ic1] = unique(aa(:,1));
        z1= zscore(C1);
        [C2,~,ic2] = unique(aa(:,2));
        z2= zscore(C2);
        for i=1:length(ic1)
            output(i).ZValue =[z1(ic1(i))  z2(ic2(i))]; 
        end
        clear aa C1 C2 ic1 ic2 z1 z2 i
        
        % find group indices in which this electrod was valid for analysis;
        % this is a sanity check! for electrodes which get to this point of
        % analysis, validGrpForThisCh must be 1, ..., 9
        validGrpForThisCh = find(elecCndNumGoodTrials(:,ch,ses) >= minElecRepPerCnd)';
        
        gr = 0;
        ch
        while gr<length(validGrpForThisCh) % for all groups in which this electrode was valid for analysis (i.e., 1,...,9)
            gr = gr+1;
            grpIdx = validGrpForThisCh(gr);
               
            % find trials within this group for which the electrode is valid for analysis
            thisChValGrpTrIdx = find( ~ismember(output(grpIdx).TrialIdx, elecCndBadTrials{gr,ch,ses}) );
            
            % get the GLM input data related to this group: mean power in the time and freq. interval of interest, separately for each trial
            thisChGLMdata{gr,ch} = ...
                      cellfun(  @(x) ...
                      squeeze( mean( reshape( Z_freq_powspctrm(x(1):x(2), x(3):x(4), output(grpIdx).TrialIdx(thisChValGrpTrIdx)), 1, [], length(thisChValGrpTrIdx) ) ) ), ...
                      mat2cell( allFreqTimIdx, ones(size(allFreqTimIdx,1),1), size(allFreqTimIdx,2) ), ...
                      'UniformOutput', 0);
            % make the predictors for this group based on number of channel valid trials in the group       
            thisChGLMpredictor{gr,ch} = repmat(output(grpIdx).ZValue,  length(thisChValGrpTrIdx) , 1);
                  
            clear grpIdx thisChValGrpTrIdx        
        end
        clear gr validGrpForThisCh freq output event Z_freq_powspctrm 
        
        % apply GLM on the current channel data
        for itvl = 1:size(allFreqTimIdx,1)
            x = cellfun(@(x) x(itvl), thisChGLMdata(:,ch));
            [glmB{ses,ch,itvl},~,glmStats{ses,ch,itvl}] = glmfit( cell2mat(thisChGLMpredictor(:,ch)), cell2mat(x), 'normal' );
            clear x
        end
        clear itvl
        
    end  % end of ch=1:...
    
    clear sesChFileNames
end
%% save the workspace
save([save_result_dir 'GLMonElecDays_workspace'])

%%  Plot histogram of the GLM coefficients across electrodes and sessions
numRows = size(allFreqTimIdx,1);
numClmn = 2;  % two arrays

figure('Name', [Monkey ', histogram of all sample coefficients'])
for rw=1:numRows
    for cl = 1:numClmn
        subplot(numRows, numClmn, ((rw-1)*numRows)+cl), hold on, box on
        title([Array_pos{cl}])
        
        % extract electrodes of Array cl
        thisArrayElecs = find(arrayIdx == cl);
        aa = glmB(:,thisArrayElecs,rw);  % extract coefficients for selected channels and freq. band
        aa = aa(:)';
        aa = cell2mat(aa);
        binEdge = -0.5:0.05:0.5;
        binCenter = mean([binEdge(1:end-1); binEdge(2:end)]);
        binCountsEv = histcounts(aa(2,:), binEdge);
        binCountsVar = histcounts(aa(3,:), binEdge);
        % plot histogram for EV coefficients
        bar(binCenter, binCountsEv/sum(binCountsEv), 'FaceColor', [153 153 255]/255);
        text(0.1, 0.2, ['\color[rgb]{' num2str([153 153 255]/255) '} expected value'])
        % plot histogram for Variance coefficients
        bar(binCenter, -binCountsVar/sum(binCountsVar), 'FaceColor', [255 153 153]/255);
        text(0.1, -0.2, ['\color[rgb]{' num2str([255 153 153]/255) '} variance '])
        xlabel('coefficient value')
        ylabel('normalized num. of samples')
    end
end

%%  Plot scatter plot of the GLM coefficients across electrodes and sessions
numRows = size(allFreqTimIdx,1);
numClmn = 2;  % two arrays

figure('Name', [Monkey ', scatter plot of the two coefficients from all samples'])
for rw=1:numRows
    for cl = 1:numClmn
        subplot(numRows, numClmn, ((rw-1)*numRows)+cl), hold on, box on
        title([Array_pos{cl}])
        
        % extract electrodes of Array cl
        thisArrayElecs = find(arrayIdx == cl);
        aa = glmB(:,thisArrayElecs,rw);  % extract coefficients for selected channels and freq. band
        aa = aa(:)';
        aa = cell2mat(aa);
        scatter(aa(2,:), aa(3,:), 'o')
        set(gca, 'XLim', [-0.25 0.25], 'YLim', [-0.25 0.25])
        xlabel('EV coefficient')
        ylabel('Var coefficient')
        clear aa thisArrayElecs
    end
end
%%  Plot histogram of the significant GLM coefficients across electrodes and sessions
numRows = size(allFreqTimIdx,1);
numClmn = 2;  % two arrays

figure('Name', [Monkey ', histogram of significant coefficients'])
for rw=1:numRows
    for cl = 1:numClmn
        subplot(numRows, numClmn, ((rw-1)*numRows)+cl), hold on, box on
        title([Array_pos{cl}])
        allSample2 = [];
        allSample3 = [];
        % extract electrodes of Array cl
        thisArrayElecs = find(arrayIdx == cl);
        % extract coefficients for selected channels and freq. band
        aa = glmB(:,thisArrayElecs,rw);  
        aa = aa(:)';
        % extract statistics for selected channels and freq. band
        sts = glmStats(:,thisArrayElecs,rw);   
        sts = sts(:)';
        % remove empty cells from the sts and aa
        tmpIdx = cellfun(@(x) isempty(x), sts, 'UniformOutput', 0); 
        sts(cell2mat(tmpIdx)) = [];
        aa(cell2mat(tmpIdx)) = [];
        clear tmpIdx
        % get the p value from statistics structure for every sample
        stsP = cellfun(@(x) getfield(x, 'p'), sts, 'UniformOutput', 0); 
        clear sts
        % find the fdr threshold for every sample
        [fdrThr,~,~] = cellfun(@(x) fdr(x), stsP, 'UniformOutput', 0);
        % find which pvalues pass the fdr threshold
        passFDR = cellfun(@(x,y) find(x<=y), stsP, fdrThr, 'UniformOutput', 0);
        for i=1:length(passFDR)
            if (~isempty(passFDR{i}) && ~isempty(find(passFDR{i} == 2)))
                allSample2 = [allSample2; aa{i}(2)];
            end
            if (~isempty(passFDR{i}) && ~isempty(find(passFDR{i} == 3)))
                allSample3 = [allSample3; aa{i}(3)];
            end
        end
        clear i aa passFDR fdrThr
        binEdge = -0.5:0.05:0.5;
        binCenter = mean([binEdge(1:end-1); binEdge(2:end)]);
        binCountsEv = histcounts(allSample2, binEdge);
        binCountsVar = histcounts(allSample3, binEdge);
        % plot histogram for EV coefficients
        bar(binCenter, binCountsEv/sum(binCountsEv), 'FaceColor', [153 153 255]/255);
        text(0.1, 0.2, ['\color[rgb]{' num2str([153 153 255]/255) '} expected value'])
        % plot histogram for Variance coefficients
        bar(binCenter, -binCountsVar/sum(binCountsVar), 'FaceColor', [255 153 153]/255);
        text(0.1, -0.2, ['\color[rgb]{' num2str([255 153 153]/255) '} variance '])
        xlabel('coefficient value')
        ylabel('normalized num. of samples')
    end
end

