clear all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder bhv_data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
bhv_data_folder             = 'Z:\data\Risk\Behavior\';
         
Monkey = 'Mojo';  % 'Moj' or 'Mac'
groupingCriteria       = 'expecter_reward & RewardVariance' ;  % $$ change the save_dir if you changed this
numConditions = 9; 
Fs = 1000;   % assume 1K sampling rate
num_ch = 96;

switch Monkey
    case 'MacDuff'
        eventTimeInterval = {[0.75 0.95]};
        alignEvent = {'cue'};  
        timeOfInt = [{[0 0.4]} , {[0 0.4]} , {[0  0.4]} , {[0.4  0.8]} , {[0.4  0.8]} , {[0.4  0.8]}];  % time in sec. based on Cue alignment
        freqOfInt = [{[4  8 ]} , {[8  14]} , {[14 30 ]} , {[4    8  ]} , {[8    14 ]} , {[14   30 ]}]; % frequency in hertz
    case 'Mojo'
        eventTimeInterval = {[0.4 0.8]};
        alignEvent = {'cue'};  
        timeOfInt = [{[0 0.4]} , {[0 0.4]} , {[0  0.4]} , {[0.4  0.8]} , {[0.4  0.8]} , {[0.4  0.8]}];  % time in sec. based on Cue alignment
        freqOfInt = [{[4  8 ]} , {[8  14]} , {[14 30 ]} , {[4    8  ]} , {[8    14 ]} , {[14   30 ]}]; % frequency in hertz
end

if strcmp(Monkey(1:3),'Mac')
    clean_data_dir = 'Z:\data\Risk\Clean Data\';
    freq_data_dir = 'H:\Risk Task\MacDuff\';  % Bahareh PC
%     freq_data_dir = 'F:\Electrode-day power spectrum\Risk Task\MacDuff\';% Saeed PC
elseif strcmp(Monkey(1:3),'Moj')
    clean_data_dir = 'Z:\data\Risk\Clean Data\';
    freq_data_dir = 'H:\Risk Task\Mojo\';  % Bahareh PC
%     freq_data_dir = 'F:\Electrode-day power spectrum\Risk Task\Mojo\';% Saeed PC
end

ApplyLogTransformToPower = 0;

% x-fold cross validation
kFold = 5;

save_dir = 'X:\Report 17\Multinomial logostic regression classifier\';

%% read our relevant lfp file names
if 0
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
% sesTrialPower1 = cell(numConditions, size(allFreqTimIdx,1)); % size= [numGroups (numFreInterest x numTimeInterest)]
% sesTrialPower2 = cell(numConditions, size(allFreqTimIdx,1)); % size= [numGroups (numFreInterest x numTimeInterest)]
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
    
    thisSesChannels1 = cell(1,size(allFreqTimIdx,1));
    thisSesChannels2 = cell(1,size(allFreqTimIdx,1));
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
        event = recursive_read(freq, 'event');  % this is the same event structure which were saved in the Clean Data folder for this session; this is going to be the same for all the channels in this sesion;
        
        % group trials based on the specified criterion
        [output] = GroupTrials(event, groupingCriteria);

        for freqtim = 1:size(allFreqTimIdx,1)
           aa = cellfun(@(x) nanmean(reshape(x,[],size(x,3)),1)'  , ...
                   cellfun(@(x) Z_freq_powspctrm( allFreqTimIdx(freqtim,1):allFreqTimIdx(freqtim,2), allFreqTimIdx(freqtim,3):allFreqTimIdx(freqtim,4), x), {output.TrialIdx}', 'UniformOutput', 0) , ...
                   'UniformOutput', 0) ;
           switch arrayIdx(ses,ch)
               case 1
                   thisSesChannels1{1,freqtim} = [thisSesChannels1{1,freqtim}, aa];
% % % %                    allSesMeanPower1(:,freqtim) = cellfun(@(x,y) [y;x] , allSesMeanPower1(:,freqtim), aa, 'UniformOutput', 0);
               case 2
                   thisSesChannels2{1,freqtim} = [thisSesChannels2{1,freqtim}, aa];
% % % %                    allSesMeanPower2(:,freqtim) = cellfun(@(x,y) [y;x] , allSesMeanPower2(:,freqtim), aa, 'UniformOutput', 0);
           end
           clear aa
        end
        
        clear freqtim output event Z_freq_powspctrm freq
    
    end  % end of ch=1:...

    save([save_dir Monkey '_Ses' num2str(ses)], 'thisSesChannels1', 'thisSesChannels2')
    clear sesChFileNames thisSesChannels1 thisSesChannels2
end

end
%%
sesFileNames = dir(save_dir);
sesFileNames = {sesFileNames.name}';
sesFileNames = sesFileNames(3:end);
sesFileNames = sesFileNames(cellfun(@(x) strcmp(x(1:3),Monkey(1:3)), sesFileNames));

grpParaVal = [1; 2; 3];
for ses = 1:length(sesFileNames) 
    display(['ses ' num2str(ses)])
    
    % load thisSesChannels1 and thisSesChannels2
    load([save_dir sesFileNames{ses}])    
    
    % convert [EV Var] grouping to Var groupign
    thisSesChannels1 = cellfun(@(x)  cellfun(@(y) cell2mat(y), mat2cell(x([1 4 7, 2 5 8, 3 6 9],:), [3 3 3], ones(1,size(thisSesChannels1{1,1},2))), 'UniformOutput', 0) , thisSesChannels1, 'UniformOutput', 0);
    thisSesChannels2 = cellfun(@(x)  cellfun(@(y) cell2mat(y), mat2cell(x([1 4 7, 2 5 8, 3 6 9],:), [3 3 3], ones(1,size(thisSesChannels2{1,1},2))), 'UniformOutput', 0) , thisSesChannels2, 'UniformOutput', 0);
    
    % get the mean power for each electrode-group
    thisSesChannels1 = cellfun(@(x) cell2mat(cellfun(@(x2) nanmean(x2), x, 'UniformOutput', 0)), thisSesChannels1, 'UniformOutput', 0);
    thisSesChannels2 = cellfun(@(x) cell2mat(cellfun(@(x2) nanmean(x2), x, 'UniformOutput', 0)), thisSesChannels2, 'UniformOutput', 0);
    
    % creat cross validation indices for each group
    xValElecIndices = crossvalind('Kfold', [1:size(thisSesChannels1{1,1},2)], kFold)';

    grpXtest = cell(1,length(thisSesChannels1));
    for k = 1:kFold
        for tim = 1:length(thisSesChannels1)
            grpXtest{1,tim}{1,k} = thisSesChannels1{1,tim}(:,find(xValElecIndices==k));  
        end
    end
    clear tim k

    piHatShfl = cell(50,kFold);
    for k = 1:kFold
        % get the test data for this round
        testData = cellfun(@(x) x(1,k), grpXtest, 'UniformOutput', 0);
        tempPara = cellfun(@(x)  cellfun(@(y) repmat(grpParaVal,1,size(x{1},2))  , x, 'UniformOutput', 0)   , testData, 'UniformOutput', 0);
        tempDataIn = cellfun(@(x) reshape(x{1,1},size(x{1,1},1)*size(x{1,1},2),1), testData, 'UniformOutput', 0);
        tempParaIn = cellfun(@(x) reshape(x{1,1},size(x{1,1},1)*size(x{1,1},2),1), tempPara, 'UniformOutput', 0);
        % remove NaN values from power samples
        testDataIn  = cellfun(@(x) x(~isnan(x)), tempDataIn, 'UniformOutput', 0);
        testDataout{k} = cellfun(@(x,y) y(~isnan(x)), tempDataIn, tempParaIn, 'UniformOutput', 0);
        clear tempPara tempPara tempDataIn tempParaIn
        
        % get the train data for this round
        aa = cellfun(@(x) x(1,find(k~=[1:kFold])), grpXtest, 'UniformOutput', 0);   
        tempPara = cellfun(@(x)  cellfun(@(y) repmat(grpParaVal,1,size(y,2))  , x, 'UniformOutput', 0)   , aa, 'UniformOutput', 0);
        
        tempDataIn  = cellfun(@(x) reshape(cell2mat(x), size(cell2mat(x),1)*size(cell2mat(x),2),1), aa, 'UniformOutput', 0);
        tempParaIn  = cellfun(@(x) reshape(cell2mat(x), size(cell2mat(x),1)*size(cell2mat(x),2),1), tempPara, 'UniformOutput', 0);
        
        % remove NaN values from power samples
        trainDataIn  = cellfun(@(x) x(~isnan(x)), tempDataIn, 'UniformOutput', 0);
        trainDataout = cellfun(@(x,y) y(~isnan(x)), tempDataIn, tempParaIn, 'UniformOutput', 0);
        clear aa tempPara tempDataIn tempParaIn
        
        % apply Multinomial logistic regression for this round
        [lfpB{k},lfpDev{k},lfpStats{k}] = cellfun(@(x,y) mnrfit(x,categorical(y),'model','nominal'), trainDataIn, trainDataout, 'UniformOutput', 0); % the model: B0 + B1*LFppower     
 
        % test accuracy of the decoder
        piHat{k} = cellfun(@(b,x) mnrval(b,x)  , lfpB{k}, testDataIn, 'UniformOutput', 0);  
        
        for shfl = 1:50
            % test accuracy of the decoder for shuffled data
            piHatShfl{shfl,k} = cellfun(@(b,x) mnrval(b,x(randperm(length(x))))  , lfpB{k}, testDataIn, 'UniformOutput', 0);  
        end
        clear shfl

% % % %         % apply Multinomial logistic regression for this round on the suffled data
% % % %         for shfl = 1:50
% % % %             [lfpShflB{shfl,k},lfpShflDev{shfl,k},lfpShflStats{shfl,k}] = cellfun(@(x,y) mnrfit(x,categorical(y(randperm(length(y)))),'model','nominal'), trainDataIn, trainDataout, 'UniformOutput', 0); % the model: B0 + B1*LFppower     
% % % %             
% % % %             % test accuracy of the decoder
% % % %             piHatShfl{shfl,k} = cellfun(@(b,x) mnrval(b,x)  , lfpShflB{shfl,k}, testDataIn, 'UniformOutput', 0);  
% % % %         end
        clear testData testDataIn trainDataIn trainDataout
    end
    clear k

    % get the decoding performance. size(decodePerf{ses}) = [k * numROIs]
    [~,predictedCategory] = cellfun(@(x)  cellfun(@(y) max(y') , x, 'UniformOutput', 0)    , piHat, 'UniformOutput', 0);
    predictedCategory = cellfun(@(x) cellfun(@(y) y', x, 'UniformOutput', 0) , predictedCategory, 'UniformOutput', 0);
    decodePerf{ses} = cell2mat(cellfun(@(x,y)  cell2mat(cellfun(@(x1,y1)  100*sum(~(x1-y1))/length(x1) , x, y, 'UniformOutput', 0))' , predictedCategory, testDataout, 'UniformOutput', 0))';

    % get the decoding performance for shuffled data. size(decodePerf{ses}) = [k * numROIs]
    [~,predictedCategoryShfl] = cellfun(@(x)  cellfun(@(y) max(y') , x, 'UniformOutput', 0)    , piHatShfl, 'UniformOutput', 0);
    predictedCategoryShfl = cellfun(@(x) cellfun(@(y) y', x, 'UniformOutput', 0) , predictedCategoryShfl, 'UniformOutput', 0);
    decodePerfShfl{ses} = cellfun(@(x,y)  cell2mat(cellfun(@(x1,y1)  100*sum(~(x1-y1))/length(x1) , x, y, 'UniformOutput', 0)) , predictedCategoryShfl, repmat(testDataout,50,1), 'UniformOutput', 0);
    decodePerfShfl{ses} = cell2mat(cellfun(@(x)  mean(cell2mat(x))'   , mat2cell(decodePerfShfl{ses}, size(decodePerfShfl{ses} ,1), ones(1,size(decodePerfShfl{ses} ,2))), 'UniformOutput', 0))';
end

    
data_noise = cell2mat(cellfun(@(xd,xn) mean(xd)-mean(xn), decodePerf, decodePerfShfl, 'UniformOutput', 0)');
errorbar(1:6, mean(data_noise), std(data_noise))
title([Monkey ', Var grouping'])
xlabel('ROIs')
ylabel('data decoding perf - noise decoding perf')
