clear
close all
clc
dbstop if error

%% Select values for thses parameters from the options otherwise the program will crash! 

Monkey            = 'Mojo'; % 'Mojo' or 'MacDuff'
groupingCriteria  = 'CuePos & TargetPos & PreGood & CurrVar (uncertainty)';
selectedArray = 'Anterior';  % 'Anterior' or 'Posterior'
% if data is already saved sapartely for arrays, set the path bellow ('if dataIsAlreadyPrep ...') where to load the data from  
dataIsAlreadyPrep = 0; % 1 if you have already saved the allSesActivityArr1 & allSesActivityArr2; otherwise 0
% define the ROIs in which you would like to plot the power. ROI is a structure with following fields
ROI.event = [{'Cue'}];
ROI.freq  = [{[10 20]}];
ROI.time  = [{[.3 .5]}];

%%
switch Monkey 
    case 'MacDuff'
        clean_data_dir  = 'Z:\data\Clean Data\';
        freq_data_dir = 'Y:\MacDuff\';
        intermediate_save_folder = 'D:\quick hack uncertainty_macduff\';
        save_dir = 'X:\Report 10\PlotGLMinput_Uncertainty_R50C50\';
        load_dir = '';
    case 'Mojo'
        clean_data_dir  = 'Z:\data\Clean Data\';
        freq_data_dir = 'Y:\Mojo\';
        intermediate_save_folder = 'D:\quick hack uncertainty_mojo\';
        save_dir = 'X:\Report 10\PlotGLMinput_Uncertainty_R50C50\';
        load_dir = '';
end

%%

Array_pos = {'Anterior','Posterior'};% Array1 is anterior and Array2 is posterior

tmp_freq_folder_list = dir(freq_data_dir);
all_freq_file_names = {tmp_freq_folder_list.name};
all_freq_file_names(1:2) = [];

tmp_clean_folder_list = dir(clean_data_dir);
all_clean_folder_names = {tmp_clean_folder_list.name};

all_clean_folder_names = all_clean_folder_names(~cellfun(@isempty, strfind(all_clean_folder_names, Monkey)));
all_clean_folder_names = all_clean_folder_names(~cellfun(@isempty, strfind(all_clean_folder_names, 'cue')));

% 16 is numOfGroups for groupingCriteria  = 'CuePos & TargetPos & PreGood & CurrVar' for C50 abd R50 conditions;
numGroups = 16;

%% load each channel-day file; zscore LFP within freq.; group (based on groupingCriteria);  and collect LFP from all sessions 
% if you have alreday collected and saved the data, do not run this section.
if ~dataIsAlreadyPrep

for ses = 1:length(all_clean_folder_names)
    tempActivity = cell(numGroups,1); 
    file_names = all_freq_file_names(~cellfun(@isempty, strfind(all_freq_file_names, all_clean_folder_names{ses}(5:end-4))));
    if isempty(file_names)
        continue
    end
    
    % sort the file names based on the channel lable index (i.e., X in the 'chX' of the freq. file name !!)
    a1 = cellfun(@(x)    str2double(x(max(strfind(x,'ch'))+2 : min(strfind(x,'.'))-1)), file_names);
    [~,sortedIdx]=sort(a1);
    file_names = file_names(sortedIdx);
    clear a1 sortedIdx
    
    
    ch = 0;
    while ch <length(file_names)
        ch = ch+1;
        % show the progress
        ['Session ' num2str(ses) ', ch ' num2str(ch)]      
        % load variable freq; size(squeeze(freq.powspctrm)) is [numTrials, numFreqBands, numTimeStamps] 
        load([freq_data_dir file_names{ch}]);
        
        % find the location of the loaded channel on the arrays
        eval(['[Array1, Array2] = electrodepinout' Monkey '(freq.label{1}(1:4));']) 
        
        % find which array this channel belongs to
        switch freq.label{1}(1:4)
            case 'chan'
                if ismember(ch, Array1)
                    arrayIdx = 1;
                else if ismember(ch, Array2)
                        arrayIdx = 2;
                    else
                        udisplay('There is an error! This channels belong to neither of the arrays!!!')
                        break
                    end
                end
                
            case 'elec'
                k = strfind(freq.label{1},'-');
                if ismember(str2num(freq.label{1}(k+1:end)), Array1)
                    arrayIdx = 1;
                else if ismember(str2num(freq.label{1}(k+1:end)), Array2)
                        arrayIdx = 2;
                    else
                        udisplay('There is an error! This channel belongs to neither of the arrays!!!')
                        break
                    end
                end
                clear k
       end
        
     if strcmp(selectedArray, Array_pos{arrayIdx})   
        % change order of freq.powspctrm dimensions and change it to [numFreqBands, numTimeStamps, numTrials]
        M0 = permute(squeeze(freq.powspctrm(:,1,:,:)), [2,3,1]);
        % turn the M0 matrix to a 2D matrix of size [numFreqBands (numTrials*numTimaStamps)]
        M = reshape(M0, size(M0,1), []);
        % zscore ezch frequency band separately 
        Z_M = (M - repmat(nanmean(M,2),1,size(M,2)))./repmat(nanstd(M,1,2),1,size(M,2));
        % reshape the zscored data into a 3D matrix of size [numFreqBands, numTimeStamps, numTrials]
        Z_freq_powspctrm = reshape(Z_M,size(M0));
      
        switch groupingCriteria
            case 'All trials'
                
            otherwise
                event = recursive_read(freq,'event');
                
                % in order to match the task conditions of the Uncertainty task with the
                % Risk task, only the following task conditions are kept:
                % C50, R50
                aa = struct2table(event);
                bb = aa.type;
                c2 = cellfun(@(x) ~isempty(strfind(x,'C50')),bb,'UniformOutput',0);
                c4 = cellfun(@(x) ~isempty(strfind(x,'R50')),bb,'UniformOutput',0);
                % find valid trila indices (for the two conditions C50 and R50)
                ind = sort([find(cell2mat(c2));find(cell2mat(c4))]);
                % make the new event structure from valid trial 
                new_event = table2struct(aa(ind,:));
                % make the new Z_freq_powspctrm from valid trials
                new_Z_freq_powspctrm = Z_freq_powspctrm(:,:,ind);
                clear aa bb c1 c2 c3 c4 ind
                
                % group remaining trials based on the groupingCriteria
                [output] = GroupTrials(new_event, groupingCriteria );% 'event' only includes successful trials, so it is not necessary to remove error trials
                for gr = 1:numel(output)
                    idx = output(gr).TrialIdx; 
                    % for the first channel in this session, initialize the variable
                    if isempty(tempActivity{gr}) 
                        tempActivity{gr} = reshape(new_Z_freq_powspctrm(:,:,idx),[1,size(new_Z_freq_powspctrm(:,:,idx))]);
                    else
                        tempActivity{gr} = cat( 1, tempActivity{gr}, reshape(new_Z_freq_powspctrm(:,:,idx),[1,size(new_Z_freq_powspctrm(:,:,idx))]) );
                    end
                end
                clear gr output event
        end  
        clear Z_freq_powspctrm new_Z_freq_powspctrm Z_M M M0
      end
    end % end in ch<...
    
    % in this session, for each trial, take the mean activity across electrodes 
    % size(tempActivity{gr}) = [#channels #fre #time #trials]
    temp = cellfun(@(x) squeeze( nanmean(x(:,:,:,:), 1) ),tempActivity, 'UniformOutput', false);     % size(temp) = [#fre #time #trials]
    save([intermediate_save_folder Monkey '_Ses' num2str(ses) '_' selectedArray],'temp')
    clear gr temp tempActivity arrayIdx

end
clear ses    

%% Integrate the intermediatly saved data and bind all the sessions  ---------------------------------------------------------
% read out file names
fNames = dir(intermediate_save_folder);
fNames(1:2) = [];
arrayIdx = find(strcmp(Array_pos, selectedArray));
outIdx = [];
% find which file names are relevantt for the selected array
for ses =1:length(fNames)
    k1 = strfind(fNames(ses).name,'_');
    k2 = strfind(fNames(ses).name,'.');
    % selected file names which should be excluded
    if ~strcmp(fNames(ses).name(k1(end)+1:k2-1), selectedArray)
        outIdx = [outIdx; ses];
    end
    clear k1 k2
end
fNames(outIdx) = [];
clear outIdx ses 

% initialize allSesActivityArr1 or allSesActivityArr2
eval(['allSesActivityArr' num2str(arrayIdx) ' = cell(numGroups,1);'])

for ses = 1:length(fNames)
    % load the session data, this will load a variable named temp which is a cell array of size 36x1
    load([intermediate_save_folder   fNames(ses).name])
    for gr=1:length(temp)
% % % %             if size(temp{gr},1)<40
% % % %                 temp{gr} = nan(47,301,2);
% % % %             end
% % % % 
        eval(['allSesActivityArr' num2str(arrayIdx) '{gr} = cat(3, allSesActivityArr' num2str(arrayIdx) '{gr}, temp{gr});'])
    end
    clear gr temp
end
clear ses

% store  allSesActivityArr1 or allSesActivityArr2 on the hard drive
save([save_dir Monkey '_allSesActivityArr' num2str(arrayIdx)], ['allSesActivityArr' num2str(arrayIdx)]) 

else
    % find the array index
    arrayIdx = find(strcmp(Array_pos,selectedArray));  % Anterior->1 ; Posterior->2
    % load the saved data for the selected array
    load([load_dir Monkey '_allSesActivityArr' num2str(arrayIdx)], ['allSesActivityArr' num2str(arrayIdx)])
end  % end of 'if ~dataIsAlreadyPrep ...'

%%
% load stored power spectrum file of one channel (any channel from any session) to regenerate the groups and their parameter values
tmp_freq_folder_list = dir(freq_data_dir);
freq_file_names = tmp_freq_folder_list(3).name;
clear tmp_freq_folder_list
load ([freq_data_dir  freq_file_names])

% get the parameter values for each group
event = recursive_read(freq,'event');
% keep only C50, R50 task conditions are kept:
aa = struct2table(event);
bb = aa.type;
c2 = cellfun(@(x) ~isempty(strfind(x,'C50')),bb,'UniformOutput',0);
c4 = cellfun(@(x) ~isempty(strfind(x,'R50')),bb,'UniformOutput',0);
% find valid trila indices (for the two conditions C50 and R50)
ind = sort([find(cell2mat(c2));find(cell2mat(c4))]);
% make the new event structure from valid trial 
new_event = table2struct(aa(ind,:));
[output] = GroupTrials(new_event, groupingCriteria);
% do not clear freq. it is needed bellow
clear event new_event aa bb c2 c4
    
% ---------------------------------------------------------------------------------------------------------------------------
roi = 0;
eval(['numFreq = size(allSesActivityArr' num2str(arrayIdx) '{1},1);'])
eval(['numTime = size(allSesActivityArr' num2str(arrayIdx) '{1},2);'])
ROI.power = cell(length(ROI.event),1);
ROI.parameters = cell(length(ROI.event),1);
while roi < length(ROI.event) % for each frequency
    roi = roi+1;
    % find time and freq interval for this ROI
    [~,freqIndex(1)] = min(abs(freq.freq - ROI.freq{roi}(1)));
    [~,freqIndex(2)] = min(abs(freq.freq - ROI.freq{roi}(2)));
    [~,timeIndex(1)] = min(abs(freq.time - ROI.time{roi}(1)));
    [~,timeIndex(2)] = min(abs(freq.time - ROI.time{roi}(2)));
    % collect data and predictor values for GLM, separately for arrays
    gr = 0;
    glmData = [];
    glmPred = [];
    while gr<length(output)
        gr = gr+1;
        eval(['glmData = [glmData; mean( squeeze(mean(allSesActivityArr' num2str(arrayIdx) '{gr}(freqIndex(1):freqIndex(2), timeIndex(1):timeIndex(2), :))) )''];'])  % some trials have NaN value which means those electrodes must be removed because of artifact
            
        eval(['glmPred = [glmPred; repmat(output(gr).Value, size(allSesActivityArr' num2str(arrayIdx) '{gr},3), 1)];'])
    end
    clear gr    
    
    ROI.power{roi} = glmData;
    ROI.parameters{roi} = glmPred;
    
    clear freqIndex timeIndex glmData glmPred
end
clear roi

%% 3D Plot the averaged power for differenc EV and Var

paraNames = {'CuePos', 'TargetPosition', 'PreGood', 'CurrVar'};
colorRange = [-0.12, 0.08];
numRows = 2;
numCulomns = 2;

% for 2D plot_Power vs Var
fg1= figure('Name', [Monkey ', ' selectedArray ', Ave. power']);
colorsString = [{'{0 1 1}'},...
                {'{0 0 1}'},...
                {'{0 0 0}'}];
colorsNumeric = [{[0 1 1]},...
                 {[0 0 1]},...
                 {[0 0 0]}];

roi = 0;
while roi < length(ROI.event) % for each frequency
    txt_2D_EV = [];
    roi = roi+1;
% % % %     figure(fg1)
% % % %     subplot(numRows, numCulomns, roi), hold on, box on
% % % %     title([ {['ROI: [' num2str(ROI.freq{roi}(1)) ' ' num2str(ROI.freq{roi}(2)) '] Hz']},{['   ' ROI.event{roi} ' [' num2str(ROI.time{roi}(1)) ' ' num2str(ROI.time{roi}(2)) '] sec']} ])
% % % %     xlabel(['Expected value'])
% % % %     ylabel(['Variance'])
% % % %     zlabel(['Ave. power in the ROI'])
% % % %     allParaValues = unique(ROI.parameters{roi}(:,3:4),'rows');  % only extract all combinations of the EV and Var
% % % %     aa = mat2cell(allParaValues,ones(size(allParaValues,1),1),2);
% % % %     % for each pair of EV and Var, extract the samples 
% % % %     [sampleIdx,~] = cellfun(@(x) ismember(ROI.parameters{roi}(:,3:4),x,'rows'), aa, 'UniformOutput', 0);
% % % %     % get the mean power for each pair of EV and Var
% % % %     powerMean = cellfun(@(x) nanmean(ROI.power{roi}(find(x))), sampleIdx, 'UniformOutput', 0);
% % % %     % get SEM of the power for each pair of EV and Var
% % % %     powerSEM = cellfun(@(x) nanstd(ROI.power{roi}(find(x)))/sqrt(sum(~isnan(ROI.power{roi}(find(x))))), sampleIdx, 'UniformOutput', 0);
    
% % % %     % plot the chart surface
% % % %     dataChart = {};
% % % %     errorChart = {};
% % % %     xData = [];
% % % %     yData = [];
% % % %     for i=1:size(allParaValues,1)
% % % %         dataChart{roi}(allParaValues(i,1)/3 , 1+sqrt(allParaValues(i,2))) = powerMean{i};
% % % %         errorChart{roi}(allParaValues(i,1)/3 , 1+sqrt(allParaValues(i,2))) = powerSEM{i};
% % % %         xData(allParaValues(i,1)/3 , 1+sqrt(allParaValues(i,2))) = allParaValues(i,1)/3 ;
% % % %         yData(allParaValues(i,1)/3 , 1+sqrt(allParaValues(i,2))) = sqrt(allParaValues(i,2));
% % % %         line([xData(allParaValues(i,1)/3 , 1+sqrt(allParaValues(i,2)))   xData(allParaValues(i,1)/3 , 1+sqrt(allParaValues(i,2)))],...
% % % %              [yData(allParaValues(i,1)/3 , 1+sqrt(allParaValues(i,2)))   yData(allParaValues(i,1)/3 , 1+sqrt(allParaValues(i,2)))],...
% % % %              dataChart{roi}(allParaValues(i,1)/3 , 1+sqrt(allParaValues(i,2)))+(errorChart{roi}(allParaValues(i,1)/3 , 1+sqrt(allParaValues(i,2)))*[-1 1]),...
% % % %              'Color','k')
% % % %     end
% % % %     clear i
% % % %     surf(xData , yData , dataChart{roi})
% % % %     plot3(xData , yData , dataChart{roi}, '.k')
% % % %     
% % % %     set(gca, 'XLim', [0 4],  'XTick', [1 2 3], 'XTickLabel', [{'3'},{'6'},{'9'}], ...
% % % %              'YLim', [-1 3], 'YTick', [0 1 2], 'YTickLabel', [{'0'},{'1'},{'4'}])
% % % %     clear aa allParaValues
    
    
    % plot 2D_Power vs EV
    uniqueVar = [0 1]';
    taskCondition = {'R50, Var = 0' , 'C50, Var > 0'};
    figure(fg1)
    subplot(numRows, numCulomns, roi), hold on, box on
    title([ {['ROI: [' num2str(ROI.freq{roi}(1)) ' ' num2str(ROI.freq{roi}(2)) '] Hz']},{['   ' ROI.event{roi} ' [' num2str(ROI.time{roi}(1)) ' ' num2str(ROI.time{roi}(2)) '] sec']} ])
    ylabel(['Ave. power in the ROI'])
    allParaValues = unique(ROI.parameters{roi}(:,4),'rows');  % only extract Var
    aa = mat2cell(allParaValues,ones(size(allParaValues,1),1),1);
    % for each Var, extract the samples 
    [sampleIdx,~] = cellfun(@(x) ismember(ROI.parameters{roi}(:,4),x,'rows'), aa, 'UniformOutput', 0);
    % get the mean power for each Var
    powerMean = cellfun(@(x) nanmean(ROI.power{roi}(find(x))), sampleIdx, 'UniformOutput', 0);
    % get SEM of the power for each Var
    powerSEM = cellfun(@(x) nanstd(ROI.power{roi}(find(x)))/sqrt(sum(~isnan(ROI.power{roi}(find(x))))), sampleIdx, 'UniformOutput', 0);

    for vr=1:length(uniqueVar)
        errorbar(1, powerMean{vr}, powerSEM{vr}, '.', 'Color', colorsNumeric{vr}, 'LineWidth', 1.5)  
        eval(['txt_2D_EV = [txt_2D_EV; {[ ''\color[rgb]'' colorsString{vr} '' '' taskCondition{vr} ]}];'])
    end
    clear vr
    set(gca, 'XLim', [0 2],  'XTick', [1], 'XTickLabel', [{'mean reward time 200ms'}])
    text(1, mean(ylim), txt_2D_EV)

end
clear roi
saveas(fg1, [save_dir Monkey '_GLMinputR50C50_' selectedArray '_' date], 'fig')
saveas(fg1, [save_dir Monkey '_GLMinputR50C50_' selectedArray '_' date], 'bmp')
close(fg1)
