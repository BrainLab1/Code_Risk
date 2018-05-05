clear
close all
clc
dbstop if error

%% Select values for thses parameters from the options otherwise the program will crash! 

Monkey            = 'MacDuff'; % 'Mojo' or 'MacDuff'
groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar & PreGood (uncertainty all conditions)';
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
        save_dir = 'X:\Report 10\PlotGLMinput_Uncertainty_AllConditions\';
        load_dir = '';
    case 'Mojo'
        clean_data_dir  = 'Z:\data\Clean Data\';
        freq_data_dir = 'Y:\Mojo\';
        intermediate_save_folder = 'D:\quick hack uncertainty_mojo\';
        save_dir = 'X:\Report 10\PlotGLMinput_Uncertainty_AllConditions\';
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
numGroups = 64;

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
                
                % Only the following task conditions are kept: C0, C25, C50, C75, C100, R25, R50, R75 
                aa = struct2table(event);
                c0 = cellfun(@(x) ~isempty(strfind(x(1:2),'C0')),aa.type,'UniformOutput',0);
                c25 = cellfun(@(x) ~isempty(strfind(x(1:3),'C25')),aa.type,'UniformOutput',0);
                c50 = cellfun(@(x) ~isempty(strfind(x(1:3),'C50')),aa.type,'UniformOutput',0);
                c75 = cellfun(@(x) ~isempty(strfind(x(1:3),'C75')),aa.type,'UniformOutput',0);
                c100 = cellfun(@(x) ~isempty(strfind(x(1:end-1),'C100')),aa.type,'UniformOutput',0);
                r25 = cellfun(@(x) ~isempty(strfind(x(1:3),'R25')),aa.type,'UniformOutput',0);
                r50 = cellfun(@(x) ~isempty(strfind(x(1:3),'R50')),aa.type,'UniformOutput',0);
                r75 = cellfun(@(x) ~isempty(strfind(x(1:3),'R75')),aa.type,'UniformOutput',0);
                % find valid trila indices (for conditions C0, C25, C50, C75, C100, R25, R50, R75)
                ind = sort([ find(cell2mat(c0)); find(cell2mat(c25)); find(cell2mat(c50)); find(cell2mat(c75)); find(cell2mat(c100)); ...
                             find(cell2mat(r25)); find(cell2mat(r50)); find(cell2mat(r75))]);
                % make the new event structure from valid trial 
                new_event = table2struct(aa(ind,:));
                % make the new Z_freq_powspctrm from valid trials
                new_Z_freq_powspctrm = Z_freq_powspctrm(:,:,ind);
                clear aa c0 c25 c50 c75 c100 r25 r50 r75 ind
                
                % group remaining trials based on the groupingCriteria
                [output] = GroupTrials(new_event, groupingCriteria );% 'event' only includes successful trials, so it is not necessary to remove error trials
                % find empty groups
                emptyGrp = cellfun(@(x) isempty(x), {output.TrialIdx})';
                output(emptyGrp) = [];
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
        eval(['aa = allSesActivityArr' num2str(arrayIdx) '{gr};']);    
        if isempty(aa)
            eval(['allSesActivityArr' num2str(arrayIdx) '{gr} = temp{gr};'])
        else
            if ~isempty(temp{gr})
                eval(['allSesActivityArr' num2str(arrayIdx) '{gr} = cat(3, allSesActivityArr' num2str(arrayIdx) '{gr}, temp{gr});'])
            end
        end
        clear aa
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
% Only the following task conditions are kept: C0, C25, C50, C75, C100, R25, R50, R75 
aa = struct2table(event);
c0 = cellfun(@(x) ~isempty(strfind(x(1:2),'C0')),aa.type,'UniformOutput',0);
c25 = cellfun(@(x) ~isempty(strfind(x(1:3),'C25')),aa.type,'UniformOutput',0);
c50 = cellfun(@(x) ~isempty(strfind(x(1:3),'C50')),aa.type,'UniformOutput',0);
c75 = cellfun(@(x) ~isempty(strfind(x(1:3),'C75')),aa.type,'UniformOutput',0);
c100 = cellfun(@(x) ~isempty(strfind(x(1:end-1),'C100')),aa.type,'UniformOutput',0);
r25 = cellfun(@(x) ~isempty(strfind(x(1:3),'R25')),aa.type,'UniformOutput',0);
r50 = cellfun(@(x) ~isempty(strfind(x(1:3),'R50')),aa.type,'UniformOutput',0);
r75 = cellfun(@(x) ~isempty(strfind(x(1:3),'R75')),aa.type,'UniformOutput',0);
% find valid trila indices (for conditions C0, C25, C50, C75, C100, R25, R50, R75)
ind = sort([ find(cell2mat(c0)); find(cell2mat(c25)); find(cell2mat(c50)); find(cell2mat(c75)); find(cell2mat(c100)); ...
             find(cell2mat(r25)); find(cell2mat(r50)); find(cell2mat(r75))]);
% make the new event structure from valid trial 
new_event = table2struct(aa(ind,:));
[output] = GroupTrials(new_event, groupingCriteria);
% do not clear freq. it is needed bellow
clear aa c0 c25 c50 c75 c100 r25 r50 r75 ind
% find empty groups
emptyGrp = cellfun(@(x) isempty(x), {output.TrialIdx})';
output(emptyGrp) = [];
clear emptyGrp
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

%% 2D Plot the averaged power for differenc EV and Var

paraNames = {'CuePos', 'TargetPosition', 'CurrEV', 'CurrVar', 'PreGood'};
colorRange = [-0.12, 0.08];
numRows = 2;
numCulomns = 2;


roi = 0;
while roi < length(ROI.event) % for each frequency
    roi = roi+1;
    allParaValues = unique(ROI.parameters{roi}(:,3:4),'rows');  % only extract all combinations of the EV and Var
    aa = mat2cell(allParaValues,ones(size(allParaValues,1),1),2);
    % for each pair of EV and Var, extract the samples 
    [sampleIdx,~] = cellfun(@(x) ismember(ROI.parameters{roi}(:,3:4),x,'rows'), aa, 'UniformOutput', 0);
    % get the mean power for each pair of EV and Var
    powerMean = cellfun(@(x) nanmean(ROI.power{roi}(find(x))), sampleIdx, 'UniformOutput', 0);
    % get SEM of the power for each pair of EV and Var
    powerSEM = cellfun(@(x) nanstd(ROI.power{roi}(find(x)))/sqrt(sum(~isnan(ROI.power{roi}(find(x))))), sampleIdx, 'UniformOutput', 0);
    
    % plot power vs EV for different variance values
    % colors for 2D plot_Power vs Var
    colorsString = [{'{0 1 1}'},...
                    {'{0 0 1}'},...
                    {'{0 0 0}'}];
    colorsNumeric = [{[0 1 1]},...
                     {[0 0 1]},...
                     {[0 0 0]}];
    fg1= figure('Name', [Monkey ', ' selectedArray ', Ave. power']);
    txt_2D = [];
    figure(fg1)
    subplot(numRows, numCulomns, roi), hold on, box on
    title([ {['ROI: [' num2str(ROI.freq{roi}(1)) ' ' num2str(ROI.freq{roi}(2)) '] Hz']},{['   ' ROI.event{roi} ' [' num2str(ROI.time{roi}(1)) ' ' num2str(ROI.time{roi}(2)) '] sec']} ])
    xlabel(['Expected value'])
    ylabel(['Ave. power in the ROI'])
    uniqueVar = unique(allParaValues(:,2));
    varianceValues = {'Var = 0' , 'Var = 0.1875' , 'Var = 0.25'};
    for vr=1:length(uniqueVar)
        tmpIdx = find(uniqueVar(vr) == allParaValues(:,2));
        errorbar(allParaValues(tmpIdx,1), [powerMean{tmpIdx}], [powerSEM{tmpIdx}], '.', 'Color', colorsNumeric{vr}, 'LineWidth', 1.5)  
        plot(allParaValues(tmpIdx,1), [powerMean{tmpIdx}], 'Color', colorsNumeric{vr}, 'LineWidth', 0.5)  
        eval(['txt_2D = [txt_2D; {[ ''\color[rgb]'' colorsString{vr} '' Var = ''  num2str(uniqueVar(vr)) ]}];'])
    end
    clear vr tmpIdx
    set(gca, 'XLim', [-0.15 1.15],  'XTick', [0 .25 .5 .75 1], 'XTickLabel', [{'0','.25','.5','.75','1'}])
    text(0.25, mean(ylim), txt_2D)
    clear txt_2D uniqueVar colorsString colorsNumeric

    % plot power vs Var for different expected values
    % colors for 2D plot_Power vs EV
    colorsString = [{'{1    0.8  1}'},...
                    {'{1    0.4  0.8}'},...
                    {'{0.8  0    0.8}'},...
                    {'{0.6  0    0.6}'},...
                    {'{0    0    0}'}];
    colorsNumeric = [{[1    0.8  1]},...
                     {[1    0.4  0.8]},...
                     {[0.8  0    0.8]},...
                     {[0.6  0    0.6]},...
                     {[0    0    0]}];
    fg2= figure('Name', [Monkey ', ' selectedArray ', Ave. power']);
    txt_2D = [];
    figure(fg2)
    subplot(numRows, numCulomns, roi), hold on, box on
    title([ {['ROI: [' num2str(ROI.freq{roi}(1)) ' ' num2str(ROI.freq{roi}(2)) '] Hz']},{['   ' ROI.event{roi} ' [' num2str(ROI.time{roi}(1)) ' ' num2str(ROI.time{roi}(2)) '] sec']} ])
    xlabel(['Variance'])
    ylabel(['Ave. power in the ROI'])
    uniqueEV = unique(allParaValues(:,1));
    expectedValues = {'EV = 0' , 'EV = 0.25' , 'EV = 0.50' , 'EV = 0.75' , 'EV = 1'};
    for ev=1:length(uniqueEV)
        tmpIdx = find(uniqueEV(ev) == allParaValues(:,1));
        errorbar(allParaValues(tmpIdx,2), [powerMean{tmpIdx}], [powerSEM{tmpIdx}], '.', 'Color', colorsNumeric{ev}, 'LineWidth', 1.5)  
        plot(allParaValues(tmpIdx,2), [powerMean{tmpIdx}], 'Color', colorsNumeric{ev}, 'LineWidth', 0.5)  
        eval(['txt_2D = [txt_2D; {[ ''\color[rgb]'' colorsString{ev} '' EV = ''  num2str(uniqueEV(ev)) ]}];'])
    end
    clear ev tmpIdx
    set(gca, 'XLim', [-0.15 0.4],  'XTick', [0 0.1875 .25], 'XTickLabel', [{'0','.1875','.25'}])
    text(0.25, mean(ylim), txt_2D)
end
clear roi
saveas(fg1, [save_dir Monkey '_GLMinput_Uncertainty_All_Pow-EV' selectedArray '_' date], 'fig')
saveas(fg1, [save_dir Monkey '_GLMinput_Uncertainty_All_Pow-EV' selectedArray '_' date], 'bmp')
close(fg1)
saveas(fg2, [save_dir Monkey '_GLMinput_Uncertainty_All_Pow-Var' selectedArray '_' date], 'fig')
saveas(fg2, [save_dir Monkey '_GLMinput_Uncertainty_All_Pow-Var' selectedArray '_' date], 'bmp')
close(fg2)
