
% - this script loads the existing LFP power spectrum file for electrod-days; zscores
%   each electrode-day data within each frequency; groups trials based on 'CuePos & TargetPos & CurrEV & CurrVar'
%   in each trial takes the mean LFP power across electrodes within the session so that for each trial,
%   we will have two activity map (one for each array); data from the two arrays are analyzed separately.
%   collects the data from all sessions and applies GLM on it. It uses zscored predictor values for GLM. 
%   then plots the GLM coefficients in the time-freq. space
% - Current settings only works for Cue aligned LFP data and groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar'
%   if you wanna change these parameters you need to adjust number of time stamps and number of groups

clear
close all
clc

%% Select values for thses parameters from the options otherwise the program will crash! 

Monkey            = 'Mojo'; % 'Mojo' or 'MacDuff'
groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar';
selectedArray = 'Anterior';  % 'Anterior' or 'Posterior'
% if data is already saved sapartely for arrays, set the path bellow ('if dataIsAlreadyPrep ...') where to load the data from  
dataIsAlreadyPrep = 1; % 1 if you have already saved the allSesActivityArr1 & allSesActivityArr2; otherwise 0

%%

switch Monkey 
    case 'MacDuff'
        bhv_data_dir  = 'Z:\data\Risk\Behavior\';
        freq_data_dir = 'D:\Power Data\';
        intermediate_save_folder = 'D:\quick hack\';
        save_dir = 'X:\Report 6\GLM on LFP power\MacDuff\';
    case 'Mojo'
        bhv_data_dir  = 'Z:\data\Risk\Behavior\';
        freq_data_dir = 'D:\Mojo-PowerData\';
        intermediate_save_folder = 'D:\Mojo_quick hack\';
        save_dir = 'X:\Report 6\GLM on LFP power\Mojo\';
end

%%

Array_pos = {'Anterior','Posterior'};% Array1 is anterior and Array2 is posterior

tmp_freq_folder_list = dir(freq_data_dir);
all_freq_file_names = {tmp_freq_folder_list.name};
all_freq_file_names(1:2) = [];

tmp_bhv_folder_list = dir(bhv_data_dir);
all_bhv_folder_names = {tmp_bhv_folder_list.name};
all_bhv_folder_names(1:2) = [];

all_bhv_folder_names = all_bhv_folder_names(~cellfun(@isempty, strfind(all_bhv_folder_names, Monkey)));          
ses_num = 0;
% here is what the numbers mean: 47 is numOfFreqBands; 36 is numOfGroups for groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar';
% 1501 is numOfTimeStamps; 48 is numOfElectrodsInEachArray  
Ses_sum1 = zeros(47*36,1501*48);
Ses_sum2 = zeros(47*36,1501*48);

% 36 is numOfGroups for groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar';
allSesActivityArr1 = cell(36,1);
allSesActivityArr2 = cell(36,1);

%% If data is already prepared for analysis, load the existing collected data from hard drive based on the previos section; 
% it runs way faster if you load and analyze the two arrays separately
% if data is already saved sapartely for arrays, set the path bellow where to load the data from  

if dataIsAlreadyPrep
    % find array index: 1 Anterior; 2 Posterior;
    arrIdx = find(strcmp(Array_pos, 'Posterior'));

    switch [Monkey ', ' selectedArray]
        case 'MacDuff, Anterior'
            load 'X:\Report 6\GLM on LFP power\MacDuff_allSesActivityArr1';
        case 'MacDuff, Posterior'
            load 'X:\Report 6\GLM on LFP power\MacDuff_allSesActivityArr2';
        case 'Mojo, Anterior'
            load 'D:\Shared\Report 5\Mojo_allSesActivityArr1';
        case 'Mojo, Posterior'
            load 'D:\Shared\Report 5\Mojo_allSesActivityArr2';
    end
end
%% load each channel-day file; zscore LFP within freq.; group (based on groupingCriteria);  and collect LFP from all sessions 
% if you have alreday collected and saved the data, do not run this section.
if ~dataIsAlreadyPrep

for ses = 1:length(all_bhv_folder_names)
    tempActivity = cell(36,1);  % 36 is numOfGroups for groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar';
    file_names = all_freq_file_names(~cellfun(@isempty, strfind(all_freq_file_names, all_bhv_folder_names{ses}(5:end-4))));
    if isempty(file_names)
        continue
    end
    ch = 0;
    while ch <length(file_names)
        ch = ch+1;
        % show the progress
        ['Session ' num2str(ses) ', ch ' num2str(ch)]      
        % load variable freq; size(squeeze(freq.powspctrm)) is [numTrials, numFreqBands, numTimeStamps] 
        load([main_folder, data_folder, freq_data_dir, file_names{ch}]);
        
        % find the location of the loaded channel on the arrays
        eval(['[Array1, Array2] = electrodepinout' Monkey '(freq.label{1}(1:4));']) 
        
        % find which array this channel belongs to
        switch freq.label{1}(1:4)
            case 'chan'
                if ismember(ch, Array1)
                    arrayIdx = 1;
                else if smember(ch, Array2)
                        arrayIdx = 2;
                    else
                        udisplay('There is an error! This channels belong to neither of the arrays!!!')
                        break
                    end
                end
                % assign electrode position on the array for current channel
                eval(['plot_indx(2,ch)   = find(' 'Array' num2str(arrayIdx) ' == ch);']);
                
            case 'elec'
                k = strfind(freq.label{1},'-');
                if ismember(str2num(freq.label{1}(k+1:end)), Array1)
                    arrayIdx = 1;
                else if ismember(str2num(freq.label{1}(k+1:end)), Array2)
                        arrayIdx = 2;
                    else
                        udisplay('There is an error! This channels belong to neither of the arrays!!!')
                        break
                    end
                end
                % assign electrode position on the array for current channel
                eval(['plot_indx(2,ch)   = find(' 'Array' num2str(arrayIdx) ' == str2num(freq.label{1}(k+1:end)));']);
                clear k
       end
        
        % assign array index to this channel
        plot_indx(1,ch) = arrayIdx;
               
        
     if strcmp(selectedArray, Array_pos{str2num(arrayIdx)})   
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
                [output] = GroupTrials(event, groupingCriteria );% 'event' only includes successful trials, so it is not necessary to remove error trials
                for gr = 1:numel(output)
                    idx = output(gr).TrialIdx; 
                    % for the first channel in this session, initialize the variable
                    if isempty(tempActivity{gr}) 
                        tempActivity{gr} = reshape(Z_freq_powspctrm(:,:,idx),[1,size(Z_freq_powspctrm(:,:,idx))]);
                    else
                        tempActivity{gr} = cat( 1, tempActivity{gr}, reshape(Z_freq_powspctrm(:,:,idx),[1,size(Z_freq_powspctrm(:,:,idx))]) );
                    end
                end
                clear gr output event
        end  
        clear Z_freq_powspctrm Z_M M M0
      end
    end % end in ch<...
    
    % in this session, for each trial, take the mean activity across electrodes 
    temp = cellfun(@(x) squeeze( nanmean(x(:,:,:,:), 1) ),tempActivity, 'UniformOutput', false);
    save([intermediate_save_folder Monkey '_Ses' num2str(ses) '_' selectedArray],'temp')
    clear gr temp tempActivity true_ch_num temp1 temp2 arr1Idx arr2Idx
    
    
end
clear ses    

% Integrate the intermediatly saved data and bind all the sessions  ---------------------------------------------------------
% read out file names
fNames = dir(intermediate_save_folder);
fNames(1:2) = [];

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
clear outIdx ses fNames

% find array index: 1 Anterior; 2 Posterior;
arrIdx = find(strcmp(Array_pos, 'Posterior'));

% initialize allSesActivityArr1 or allSesActivityArr2
eval(['allSesActivityArr' num2str(arrIdx) ' = cell(36,1);'])

for ses = 1:length(fNames)
    % load the session data, this will load a variable named temp which is a cell array of size 36x1
    load([intermediate_save_folder   fNames(ses).name])
    for gr=1:length(temp)
        eval(['allSesActivityArr' num2str(arrIdx) '{gr} = cat(3, allSesActivityArr' num2str(arrIdx) '{gr}, temp{gr});'])
    end
    clear gr temp
end
clear ses

% store  allSesActivityArr1 or allSesActivityArr2 on the hard drive
save([save_dir Monkey '_allSesActivityArr' num2str(arrIdx)], ['allSesActivityArr' num2str(arrIdx)]) 

end  % end of 'if ~dataIsAlreadyPrep ...'

%% Apply GLM; it runs way faster if you load and analyze the two arrays separately
% load stored power spectrum file of one channel (any channel from any session) to regenerate the groups and their parameter values
tmp_freq_folder_list = dir(freq_data_dir);
freq_file_names = tmp_freq_folder_list(3).name;
clear tmp_freq_folder_list
load ([freq_data_dir  freq_file_names])

% get the parameter values for each group
event = recursive_read(freq,'event');
[output] = GroupTrials(event, groupingCriteria);
% do not clear freq. it is needed bellow
clear event

% zscore regressors
cuePosZ   =  [-1           1
              -0.7071   0.7071];
trgetPosZ =  [-1           1
              -0.7071   0.7071]; 
evZ       =  [3   6  9
              -1  0  1];
varZ      =  [0        1      4
              -0.80  -0.32  1.12];
for gr = 1:length(output)
    idx1 = find( output(gr).Value(1) == cuePosZ(1,:) );
    idx2 = find( output(gr).Value(2) == trgetPosZ(1,:) );
    idx3 = find( output(gr).Value(3) == evZ(1,:) );
    idx4 = find( output(gr).Value(4) == varZ(1,:) );
    output(gr).Value = [cuePosZ(2,idx1), trgetPosZ(2,idx2), evZ(2,idx3), varZ(2,idx4)];
    clear idx1 idx2 idx3 idx4
end
clear gr
    
% ---------------------------------------------------------------------------------------------------------------------------
fq = 0;
b=cell(1);
stats=cell(1); 
eval(['numFreq = size(allSesActivityArr' num2str(arrIdx) '{1},1);'])
eval(['numTime = size(allSesActivityArr' num2str(arrIdx) '{1},2);'])
while fq < numFreq % for each frequency
    fq = fq+1;
    t = 0;
    while t < numTime % for each time sample
        t = t+1;
        % collect data and predictor values for GLM, separately for arrays
        gr = 0;
        glmData = [];
        glmPred = [];
        while gr<length(output)
            gr = gr+1;
            eval(['glmData = [glmData; squeeze(allSesActivityArr' num2str(arrIdx) '{gr}(fq,t,:))];'])  % some trials have NaN value which means those electrodes must be removed because of artifact
            
            eval(['glmPred = [glmPred; repmat(output(gr).Value, length(squeeze(allSesActivityArr' num2str(arrIdx) '{gr}(fq,t,:))), 1)];'])
        end
        clear gr    

%         [b{fq,t},~,stats{fq,t}] = glmfit(glmPred, glmData-min(glmData)+eps, 'gamma');
        [b{fq,t},~,stats{fq,t}] = glmfit(glmPred, glmData, 'normal');
        [fq t]
    end
end


%% Plot the GLM coefficient for each time and frequency; use False Discovery Rate for significance level

paraNames = {'CuePos', 'TargetPosition', 'CurrEV', 'CurrVar'};
fg1= figure('Name',Array_pos{1}); fg2= figure('Name',Array_pos{2});
colorRange = [-0.12, 0.08];

for para = 2:5 % for each parameter except for the constant value!
    
    temp = cellfun(@(x) x(para), b, 'UniformOutput', 0);
    % get the p value for this parameter
    pVal = cellfun(@(x) x.p(para), stats, 'UniformOutput', 0); 
    pThr = cellfun(@(x) fdr(x.p), stats, 'UniformOutput', 0); % find the significant pThreshold defined by False Discovery Rate method
    signifP = cell2mat(pVal) <= cell2mat(pThr);
    figure(fg1)
    subplot(2,2,para-1), hold on, box on
    title(paraNames{para-1})
    H = mesh(cell2mat(temp));
    colormap('jet'),view([0 90]),grid on,axis('tight')
    caxis([min(cell2mat(temp(:)))-0.01 max(cell2mat(temp(:)))+0.01]);
    H.XData = freq.time;
    set(gca,'XTick', [0 0.4 0.8], 'XTickLabel', {'CueOn','CueOff','TargetOn'})
    line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    line([0.4,0.4],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    line([0.8,0.8],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    ylabel('Frequency')
%     [min(cell2mat(temp1(:))) max(cell2mat(temp1(:)))]
    colorbar
    contour(H.XData, H.YData, signifP, 1, 'LineColor', 'k', 'ContourZLevel', 1.5);
    
    clear temp pVal pThr signifP
end
