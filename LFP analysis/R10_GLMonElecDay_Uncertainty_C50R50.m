
% last update: 11.04.2017, Bahareh: FDR is applied for slecting significant p values. event-time outliers were removed
% last update 25.10.2017: this is checked if arrayIdx is changed during
% each session. if there is any change, an error occures. so, it is
% neccessary to redraw the results of Report 7.
% This script loads LFP and event data from the Clean Data folder; groups trials based on 'expecter_reward & RewardVariance';  
% For each electrode in each condition, determines in which trials the electrode had bad signal and has to be removed
% from the analysis; it removes all electrodes from the analysis for trials with wierd (outlier!) event timing;
% In each session, we applied GLM on those channels which had enough repetition-per-group for all groups 
% GLM input data was mean power in the time and freq. interval of interest, measured separately in each trial
% We zscored the parameter values which were used for grouping trials and added a new field .ZValue to the output. These
% values will serve as GLM predictors
 


% clear
% close all
% clc

%% set some parameters
Monkey            = 'MacDuff';
groupingCriteria  = 'type';
num_ch = 96;
minElecRepPerCnd = 10;
timeOfInt = [{[0.3 0.5]}];  % time in sec. based on Cue alignment
freqOfInt = [{[10 20]}]; % frequency in hertz

%%
if strcmp(Monkey,'MacDuff')
    clean_data_dir = 'Z:\data\Clean Data\';
%     freq_data_dir = 'H:\Risk Task\MacDuff\'; %Bahareh PC
    freq_data_dir = 'E:\Electrode-day power spectrum\MacDuff\';% Saeed PC
    save_result_dir = 'D:\Shared\Report 10\GLM on elec-days_Uncertainty_R50C50\MacDuff\';
elseif strcmp(Monkey,'Mojo')
    clean_data_dir = 'Z:\data\Clean Data\';
%     freq_data_dir = 'H:\Risk Task\Mojo\';% Bahareh PC
    freq_data_dir = 'E:\Electrode-day power spectrum\Mojo\';% Saeed PC
    save_result_dir = 'D:\Shared\Report 10\GLM on elec-days_Uncertainty_R50C50\Mojo\';
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
% for ses = 1:length(monkeyCleanData) 
%     ses
%     % load the current session data; this will load a variable named data which include both LFP signals and event info 
%     load([clean_data_dir monkeyCleanData{ses}])
%     
%     % read out the event structure
%     event = recursive_read(data,'event');
%     
%     % find trial indices for conditions "F9" and "B9H" in the event structure 
%     trlIdx2Conditions = find(cellfun(@(x) strcmp(x,'C50a')+strcmp(x,'C50b')+strcmp(x,'R50a')+strcmp(x,'R50b'),{event.type}'));
%     % keep only conditions "F9" and "B9H" from the event structure
%     event = event(trlIdx2Conditions);    
%     
%     event = struct2table(event);
%     event.type = cellfun(@(x) x(1:3),event.type,'UniformOutput',0);
%     event = table2struct(event);
%     
%     % group trials based on the specified criterion
%     [output] = GroupTrials(event, groupingCriteria);
%     
%     % convet event structure to table
%     event = struct2table(event);
%     
%     % initialize the variable if doesn't exist !! this variable tells us in each condition, for which trials the electrode  
%     % must be removed; it contains trial indices
%     if ~exist('elecCndBadTrials','var') 
%         elecCndBadTrials = cell(size(output,1), 96, length(monkeyCleanData)); % size: [numGroups x numElectrodes x numSessions]
%     end
%     
%     % this variable shows number of valid trials per condition  for each electrode
%     if ~exist('elecCndNumGoodTrials','var') 
%         elecCndNumGoodTrials = nan(size(output,1), 96, length(monkeyCleanData) ); % size: [numGroups x numElectrodes x numSessions]
%     end
%     
%     % fill in the elecCndBadTrials cells with trial numbers
%     for gr=1:size(output,1)
%         for str = 1:length(output(gr).TrialIdx) % for each trials within this group
%             trIdx = output(gr).TrialIdx(str);
%             
%             % pool bad_lfp_channels and bad_amplitude_channels
%             allBadCh = unique( [event.bad_lfp_channels{trIdx,1}'] );
%             
%             % add this trial index to for the bad channels in this trial!!!
%             elecCndBadTrials(gr,allBadCh,ses) =  cellfun( @(x) {[elecCndBadTrials{gr,x,ses}; trIdx]}, mat2cell(allBadCh,repmat(1,length(allBadCh),1),1) );
%             clear trIdx allBadCh
%         end
%         % find number of good trials in this condition for each electrode 
%         elecCndNumGoodTrials(gr,:,ses) = length(output(gr).TrialIdx) - cellfun(@length, elecCndBadTrials(gr,:,ses));
%     end
%     clear gr data event output
%    
% end

%% 
freqDataFiles = dir(freq_data_dir);
freqDataFiles(1:2) = [];
freqDataFiles = {freqDataFiles.name}';

% initialize the variable!!
allSesMeanPower = cell(length(monkeyCleanData), 9, length(freqOfInt), length(timeOfInt)); % size= [numSesssions x numGroups x numFreInterest x numTimeInterest]
glmB = cell(length(monkeyCleanData) ,96, 2);
glmStats = cell(length(monkeyCleanData) ,96, 2);

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
        
        % load the electrode power spectrum if it has passed the above criterion; this will load the freq variable ... 
        load([freq_data_dir sesChFileNames{ch}])
        
        
        % find which array this channel belongs to
        eval(['[Array1, Array2] = electrodepinout' Monkey '(freq.label{1}(1:4));']) 
        switch freq.label{1}(1:4)
            case 'chan'
                if ismember(ch, Array1)
                    arrayIdx(ses,ch) = 1;
                else if ismember(ch, Array2)
                        arrayIdx(ses,ch) = 2;
                    else
                        display('There is an error! This channels belong to neither of the arrays!!!')
                        break
                    end
                end
                
            case 'elec'
                k = strfind(freq.label{1},'-');
                if ismember(str2num(freq.label{1}(k+1:end)), Array1)
                    arrayIdx(ses,ch) = 1;
                else if ismember(str2num(freq.label{1}(k+1:end)), Array2)
                        arrayIdx(ses,ch) = 2;
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
% % % %             % get all combinations of freq and time of interest intervals -> allFreTimInterest(i,:) = [freInteresStart freInteresEnd timInterestStart timInterestEnd]
% % % %             allFreqTimIdx = combvec(fIdx',tIdx')';
            
            % pair time and frequency indices
            allFreqTimIdx = [fIdx , tIdx];
            
            allROIinfo = {};
            for tf=1:size(allFreqTimIdx,1)
                allROIinfo{tf} = ['freq[' num2str(round(freq.freq(allFreqTimIdx(tf,1)))) ',' num2str(round(freq.freq(allFreqTimIdx(tf,2)))), ...
                                  '] Cue[' num2str(freq.time(allFreqTimIdx(tf,3))) ',' num2str(freq.time(allFreqTimIdx(tf,4))) ']' ];
            
                
            end
            clear tIdx fIdx tf
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
        
        % find trial indices for conditions "F9" and "B9H" in the event structure 
        trlIdx2Conditions = find(cellfun(@(x) strcmp(x,'C50a')+strcmp(x,'C50b')+strcmp(x,'R50a')+strcmp(x,'R50b'),{event.type}'));
        % keep only conditions "F9" and "B9H" from the event structure
        event = event(trlIdx2Conditions);    
        
        event = struct2table(event);
        event.type = cellfun(@(x) x(1:3),event.type,'UniformOutput',0);
        event = table2struct(event);        
        
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

%% display an error if the array index do not have the same indices in all the sessions
if size(unique(arrayIdx,'rows'),1) >1
    error('ERROR: channel indices are different across sessions')
end

%% save the workspace
save([save_result_dir 'GLMonElecDays_workspace_' Monkey '_' date])

%%  Plot histogram of the significant GLM coefficients across electrodes and sessions (without using FDR)
numRows = size(allFreqTimIdx,1);
numClmn = 2;  % two arrays
pThr = 0.05;

figure('Name', [Monkey ', histogram of coefficients'])
for rw = 1:numRows
    for cl = 1:numClmn
        subplot(numRows, numClmn, ((rw-1)*numClmn)+cl), hold on, box on
        title({[Array_pos{cl} ', Significant: pGLM <= 0.05'], ...
                allROIinfo{rw}})
        allSampleEV = [];
        allSampleVar = [];
        % extract electrodes of Array cl
        thisArrayElecs = find(arrayIdx(1,:) == cl);
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
        stsP = cell2mat(stsP);
        % find which coefficients are individually significant (GLM pvalue<pThr)
        significantCoeff = stsP <= pThr;
        % get the coefficient values
        aa = cell2mat(aa);
        allSampleEV = aa(2,:);
        allSampleVar = aa(3,:);
        clear aa 
        % define bins!
        binEdge = -0.5:0.01:0.5;
        binCenter = mean([binEdge(1:end-1); binEdge(2:end)]);
        % get the histogram of all samples
        binCountsEV = histcounts(allSampleEV, binEdge);
        binCountsVar = histcounts(allSampleVar, binEdge);
        % get the histogram of significant samples
        signiBinCountsEv = histcounts(allSampleEV(significantCoeff(2,:)), binEdge);
        signiBinCountsVar = histcounts(allSampleVar(significantCoeff(3,:)), binEdge);
        % plot histogram for all EV and Var coefficients
        bar(binCenter, binCountsEV/sum(binCountsEV), 'FaceColor', [170 170 170]/255, 'EdgeColor', 'none');
        bar(binCenter, -binCountsVar/sum(binCountsVar), 'FaceColor', [170 170 170]/255, 'EdgeColor', 'none');
        % plot histogram for significant EV and Var coefficients
        bar(binCenter, signiBinCountsEv/sum(binCountsEV), 'FaceColor','b', 'EdgeColor', 'none');
        bar(binCenter, -signiBinCountsVar/sum(binCountsVar), 'FaceColor', 'm', 'EdgeColor', 'none');
        
        text(0.05, 0.15,  ['\color[rgb]{'  num2str([0 0 1]) '} expected value'])
        text(0.05, 0.1, ['\color[rgb]{' num2str([0 0 1]) '} significants ' num2str(sum(signiBinCountsEv)) '/' num2str(sum(binCountsEV))])

        text(0.05, -0.1,  ['\color[rgb]{'  num2str([1 0 1]) '} variance'])
        text(0.05, -0.15, ['\color[rgb]{' num2str([1 0 1]) '} significants ' num2str(sum(signiBinCountsVar)) '/' num2str(sum(binCountsVar))])

        % add more info to the plot
        set(gca, 'YLim', [-0.2 0.2], 'YTick', [-0.2 -0.1 0 0.1 0.2], 'YTickLabel', {'20', '10', '0', '10', '20'}, 'XTick', [-0.4 :0.2: 0.4], 'XTickLabel', {'-0.4', '-0.2', '0', '0.2', '0.4'})
        plot([0 0], [-.2 .2], 'k--')
        plot([-0.5 0.5], [0 0], 'k')
        xlabel('coefficient value')
        ylabel('percentage of samples')
    end
end

