
% this script applies GLM on LFP (in pre-defined frequency and time interval); applies GLM on the lick signal (in pre-defined
% time intervals, could be more than one interval); applies GLM on Z(Log(RT)); then correlates GLM coefficients for EV and Var
% for LFP-Lick and LFP-RT; both LFP and event-time outliers are removed

clear
clc
dbstop if error

%%

Monkey           = 'Mojo';% 'Mojo' or 'MacDuff'
groupingCriteria = 'expecter_reward & RewardVariance';
minElecRepPerCnd = 1;
alignEvent       = [   {'cue'}  ,  {'reward'} ]; % select between 'cue' and 'reward'
timeOfIntLick    = [{[0.4 0.8]} , {[-.377 0]} ];
timeOfIntSpctrm  = [{[0.2 0.5]}];  % time in sec. based on Cue alignment
freqOfIntSpctrm  = [{[10 20]}]; % frequency in hertz
Fs               = 1000;
clean_data_dir  = 'Z:\data\Risk\Clean Data\';
% freq_data_dir   = ['F:\Electrode-day power spectrum\' Monkey '\'];
freq_data_dir   = ['H:\' Monkey '\'];
save_result_dir = ['X:\Report 8\LickLFPCorr\' Monkey '\'];

Array_pos = {'Anterior','Posterior'};% Array1 is anterior and Array2 is posterior

% read out all file names in the clean data folder
cleanDataFiles = dir(clean_data_dir);
cleanDataFiles(1:2) = [];
cleanDataFiles = {cleanDataFiles.name}';

% select file names related to the selected monkey
monkeyCleanData = cleanDataFiles(~cellfun(@isempty, strfind(cleanDataFiles, Monkey)));  

%% Apply GLM on power, lick signal and reactiont time (session-by-session); exclude all the outliers

freqDataFiles = dir(freq_data_dir);
freqDataFiles(1:2) = [];
freqDataFiles = {freqDataFiles.name}';

sesLickSig1 = {[]};
sesLickSig2 = {[]};

sesTrialChPowROI(1,1:96,1:length(monkeyCleanData)) = {NaN};

for ses = 1:length(monkeyCleanData) 
    ses

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
    
    % initialize arrayIdxfor this sesion
    arrayIdx = nan(1,96);
    
    for ch = 1:length(sesChFileNames)
        [ses ch]
        % load the electrode power spectrum if it has passed the above criterion; this will load the freq variable ... 
        load([freq_data_dir sesChFileNames{ch}])
        
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
        
        
        % find indices for frequency and time of interest
        if ~exist('allFreqTimIdx','var')
            % find the indices for desired freq 
            [~,a] = cellfun(@(x) min(abs(x(1)-freq.freq)) , freqOfIntSpctrm);
            fIdx(:,1) = a';
            clear a
            [~,a] = cellfun(@(x) min(abs(x(2)-freq.freq)) , freqOfIntSpctrm);
            fIdx(:,2) = a';
            clear a
            % find the indices for desired time 
            [~,a] = cellfun(@(x) min(abs(x(1)-freq.time)) , timeOfIntSpctrm);
            tIdx(:,1) = a';
            clear a
            [~,a] = cellfun(@(x) min(abs(x(2)-freq.time)) , timeOfIntSpctrm);
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
        
        for tr = 1:length(event)            

        sesTrialChPowROI{tr,ch,ses} = ...
            cell2mat(...
                cellfun(  @(x) ...
                squeeze( mean( reshape( Z_freq_powspctrm(x(1):x(2), x(3):x(4), tr), 1, [], 1 ) ) ), ...
                mat2cell( allFreqTimIdx, ones(size(allFreqTimIdx,1),1), size(allFreqTimIdx,2) ), ...
                'UniformOutput', 0)...
                    );
                  
        end
        clear Z_freq_powspctrm tr event freq
        
    end % channel counter end
    clear ch
% % % %     
    
    % initialize sesLickGLM and sesLFPGLM and sesRTGLM variables for this sesion
    for lickIntv = 1:length(timeOfIntLick)
        sesLickGLM{ses, lickIntv} = struct('bArray1',[],'pArray1',[],'bArray2',[],'pArray2',[]); % the cell includes a structure with filds barray1 bArray2 pArray1 pArray2
    end
    clear lickIntv
    sesLFPGLM{ses,1} = struct('bArray1',[],'pArray1',[],'bArray2',[],'pArray2',[]); % the cell includes a structure with filds barray1 bArray2 pArray1 pArray2
    sesRTGLM{ses,1} = struct('bArray1',[],'pArray1',[],'bArray2',[],'pArray2',[]); % the cell includes a structure with filds barray1 bArray2 pArray1 pArray2

    for ary = 1:2
        % find electrodes on this array
        a = sesTrialChPowROI(:,find(arrayIdx == ary),ses);
        % for each trial get the averaged power on the array; remove the NaN values from the invalid lectrodes
        for tr = 1:size(sesTrialChPowROI(:,:,ses),1) % for each trial
            meanPowTrlArray{tr,1} = nanmean(cell2mat(a(tr,:)),2);
        end
        clear tr a
        meanPowTrlArray = cell2mat(meanPowTrlArray);
        % indeces of trials with NaN values in power spectrum are stored in 'nanInd'; These are trials for which all electrodes were invalid 
        nanInd = find(isnan(meanPowTrlArray(:,1)));
        
        % load channel 96 of this session to get the session event table
        load([freq_data_dir sesChFileNames{96}])
        % read out the event structure
        event = recursive_read(freq, 'event');  % this is the same event structure which were saved in the Clean Data folder for this session
        
        % ---- LFP outlier removal: remove trials for which there is no valid electrode in this array from both the event structure and collective power
        meanPowTrlArray(nanInd,:) = [];
        event(nanInd) = [];
        % ---- Event timing & RT outlier removal: remove trials for which event times or reaction times were outliers
        evtTimeAndRTOutTrials = cellfun(@(x,y,z) (x==1)||(y==1)||(z==1), {event.OutCueOnOffZ2}', {event.OutAqToRwdZ2}', {event.OutRTZ2}');
        meanPowTrlArray(evtTimeAndRTOutTrials,:) = [];
        event(evtTimeAndRTOutTrials) = [];
        
        % get the mean lick signal in the desired time interval for the remaining trials
        for tr = 1:length(event) % for each trial
            for lickIntv = 1:length(timeOfIntLick)
                switch alignEvent{lickIntv}
                    case 'cue'
                        cutInterval = event(tr).DiodeEventTime.CueOnset + (timeOfIntLick{lickIntv}*Fs);
                        tmpSig = event(tr).LickSignal > 2.5;
                        temp = tmpSig(cutInterval(1) : cutInterval(2))';
                        clear tmpSig cutInterval
                    case 'reward'
                        cutInterval = event(tr).RewardOnTime(1) + (timeOfIntLick{lickIntv}*Fs);
                        tmpSig = event(tr).LickSignal > 2.5;
                        temp = [tmpSig(cutInterval(1) : cutInterval(2))'];
                        clear tmpSig cutInterval
                end
                eval(['sesLickSig' num2str(ary) '{ses}(tr,lickIntv) = mean(temp);'])
                clear temp
            end
            
        end
        clear tr lickIntv
        
        % get the reaction time [i.e., zscore(log(RT))] for the remaining trials
        eval(['sesRT' num2str(ary) '{ses,1}= cellfun(@(x) x, {event.z_log_RT})'';'])
        
        % for the current array and session, group trials based on the specified criterion to have the grouping parameter values
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
    
        % for this session, initialize the GLM predictor matrix with NaN
        eval(['sesGLMPred' num2str(ary) '{ses} = nan(length(event),size(output(1).ZValue, 2));'])
        % for this session, fill the predictor matrix using output.ZValue
        for gr = 1:length(output)
            eval(['sesGLMPred' num2str(ary) '{ses}(output(gr).TrialIdx,:) = repmat(output(gr).ZValue, length(output(gr).TrialIdx), 1);'])
        end     
        clear gr 
        
        % apply GLM on lick signal for the valid trials of this array
        for lickIntv = 1:length(timeOfIntLick)
            eval(['[b,~,glmStats] = glmfit( sesGLMPred' num2str(ary) '{ses}, sesLickSig' num2str(ary) '{ses}(:,lickIntv), ''normal'' );'])
            eval(['sesLickGLM{ses,lickIntv}.bArray' num2str(ary) ' = b;'])
            eval(['sesLickGLM{ses,lickIntv}.pArray' num2str(ary) ' = glmStats.p;'])
            clear glmStats b
        end
        
        % apply GLM on LFP power for the valid trials of this array
        eval(['[b,~,glmStats] = glmfit( sesGLMPred' num2str(ary) '{ses}, meanPowTrlArray, ''normal'' );'])
        eval(['sesLFPGLM{ses}.bArray' num2str(ary) ' = b;'])
        eval(['sesLFPGLM{ses}.pArray' num2str(ary) ' = glmStats.p;'])
        clear glmStats b

        % apply GLM on Reaction Times 
        eval(['[b,~,glmStats] = glmfit( sesGLMPred' num2str(ary) '{ses}, sesRT' num2str(ary) '{ses,1}, ''normal'' );'])
        eval(['sesRTGLM{ses,1}.bArray' num2str(ary) ' = b;'])
        eval(['sesRTGLM{ses,1}.pArray' num2str(ary) ' = glmStats.p;'])
        clear glmStats b
        
        clear glmStats meanPowTrlArray  event output
        
    end % end of ary = 1:2...
    clear ary    
        
    clear sesChFileNames
end
clear ses    

save([Monkey '_Workspace'])
%% Plot the results


for lickIntv = 1:length(timeOfIntLick)
    fg = figure('Name', [Monkey ', Lick signal in ' alignEvent{lickIntv} ' [' num2str(timeOfIntLick{lickIntv}(1)) ' ' num2str(timeOfIntLick{lickIntv}(2)) ']']); 
    numRows =2; % this is number of GLM factors
    glmFactorName = [{'EV'},{'Var'}];
    numColumns = 2; % this is number of arrays
    
    for rw = 1:numRows
        for clm = 1:numColumns
            subplot(numRows, numColumns, ((rw-1)*numColumns)+clm), hold on, box on
            title([Array_pos{clm} ', ' glmFactorName{rw} ', Lick interval = ' alignEvent{lickIntv} ' [' num2str(timeOfIntLick{lickIntv}(1)) ',' num2str(timeOfIntLick{lickIntv}(2)) ']'])
            xlabel('Lick GLM coefficient')
            ylabel('LFP GLM coefficient')
            eval(['cellfun(@(x,y) scatter(x.bArray' num2str(clm) '(rw+1), y.bArray' num2str(clm) '(rw+1), ''b'', ''filled''), sesLickGLM(:,lickIntv) , sesLFPGLM )'])
            plot([-0.2 0.2], [-.2 .2], '-', 'Color', [128 128 128]/255)
            eval(['x = cellfun(@(x,y) [x.bArray' num2str(clm) '(rw+1), y.bArray' num2str(clm) '(rw+1)], sesLickGLM(:,lickIntv) , sesLFPGLM, ''UniformOutput'', 0);'])
            x = cell2mat(x);
% % % %             [regSlope,regSlopInt,~,~,regStats] = regress(x(:,2), [x(:,1) ones(length(x(:,2)),1)]);
% % % %             xfit = -.2:0.05:.2;
% % % %             ymodel = regSlope(1)*xfit + regSlope(2);
% % % %             plot(xfit, ymodel, '--r')
% % % %             text(-0.18, 0.1, [{'fitting method : MATLAB regress'} , {['Slope = ' num2str(regSlope(1))]} , {['Slope interval = [' num2str(regSlopInt(1,1)) ' , ' num2str(regSlopInt(1,2)) ']']} , {[' R^2 = ' num2str(regStats(1))]}] )
% % % %             clear x regSlope regSlopInt regStats xfit ymodel
            [robustSlope,robustStats] = robustfit(x(:,1), x(:,2));
            xfit = -.2:0.05:.2;
            ymodel = robustSlope(1) + robustSlope(2)*xfit;
            plot(xfit, ymodel, '--r')
            text(-0.18, 0.1, [{'fitting method : MATLAB robustfit'} , {['Slope = ' num2str(robustSlope(2))]} , {['Slope interval = [' num2str(robustSlope(2)-robustStats.se(2)) ' , ' num2str(robustSlope(2)+robustStats.se(2)) ']']}] )
            clear x  xfit ymodel
            % devide the 4 quadrants
            plot([-.2 .2],[0 0], [0 0],[-.2 .2], 'Color', [128 128 128]/255)
        end
    end

end
    
    

fg = figure('Name', [Monkey ', Reaction Time and LFP']); 
numRows =2; % this is number of GLM factors
glmFactorName = [{'EV'},{'Var'}];
numColumns = 2; % this is number of arrays

for rw = 1:numRows
    for clm = 1:numColumns
        subplot(numRows, numColumns, ((rw-1)*numColumns)+clm), hold on, box on
        title([Array_pos{clm} ', ' glmFactorName{rw}])
        xlabel('Reaction Time GLM coefficient')
        ylabel('LFP GLM coefficient')
        eval(['cellfun(@(x,y) scatter(x.bArray' num2str(clm) '(rw+1), y.bArray' num2str(clm) '(rw+1), ''b'', ''filled''), sesRTGLM , sesLFPGLM )'])
        plot([-0.2 0.2], [-.2 .2], '-', 'Color', [128 128 128]/255)
        eval(['x = cellfun(@(x,y) [x.bArray' num2str(clm) '(rw+1), y.bArray' num2str(clm) '(rw+1)], sesRTGLM , sesLFPGLM, ''UniformOutput'', 0);'])
        x = cell2mat(x);
        [robustSlope,robustStats] = robustfit(x(:,1), x(:,2));
        xfit = -.2:0.05:.2;
        ymodel = robustSlope(1) + robustSlope(2)*xfit;
        plot(xfit, ymodel, '--r')
        text(-0.18, 0.1, [{'fitting method : MATLAB robustfit'} , {['Slope = ' num2str(robustSlope(2))]} , {['Slope interval = [' num2str(robustSlope(2)-robustStats.se(2)) ' , ' num2str(robustSlope(2)+robustStats.se(2)) ']']}] )
        clear x  xfit ymodel
        % devide the 4 quadrants
        plot([-.2 .2],[0 0], [0 0],[-.2 .2], 'Color', [128 128 128]/255)
    end
end


