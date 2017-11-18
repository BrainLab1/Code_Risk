
% - this script loads the existing LFP power spectrum file for electrod-days; zscores
%   each electrode-day data within each frequency; groups trials based on 'CuePos & TargetPos & CurrEV & CurrVar'
%   in each trial takes the mean LFP power across electrodes within the session so that for each trial,
%   we will have two activity map (one for each array); for each session trials with outlier event times are removed
%   when the event data is read out for each channel, outlier trials are removed from channel power and event. 
%   data from the two arrays are analyzed separately.
%   collects the data from all sessions and applies GLM on it. It uses zscored predictor values for GLM. 
%   then plots the GLM coefficients in the time-freq. space
% - Current settings only works for Cue aligned LFP data and groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar'
%   if you wanna change these parameters you need to adjust number of time stamps and number of groups

clear
close all
clc

%% Select values for thses parameters from the options otherwise the program will crash! 

Monkey            = 'Mojo'; % 'Mojo' or 'MacDuff'
groupingCriteria  = 'CuePos & TargetPos & PreGood & CurrEV (uncertainty)';
selectedArray = 'Anterior';  % 'Anterior' or 'Posterior'
% if data is already saved sapartely for arrays, set the path bellow ('if dataIsAlreadyPrep ...') where to load the data from  
dataIsAlreadyPrep = 0; % 1 if you have already saved the allSesActivityArr1 & allSesActivityArr2; otherwise 0

%%

switch Monkey 
    case 'MacDuff'
        clean_data_dir  = 'Z:\data\Clean Data\';
        freq_data_dir = 'E:\Electrode-day power spectrum\MacDuff';
        intermediate_save_folder = 'D:\quick hack_macduff\';
        save_dir = 'D:\Shared\Report 9\GLM analysis\MacDuff\GLM input data\';
    case 'Mojo'
        clean_data_dir  = 'Z:\data\Clean Data\';
        freq_data_dir = 'E:\Electrode-day power spectrum\Mojo\';
        intermediate_save_folder = 'D:\quick hack_Mojo\';
        save_dir = 'D:\Shared\Report 9\';
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

ses_num = 0;
% here is what the numbers mean: 47 is numOfFreqBands; 40 is numOfGroups for groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar';
% 301 is numOfTimeStamps; 48 is numOfElectrodsInEachArray  
Ses_sum1 = zeros(47*40,301*48);
Ses_sum2 = zeros(47*40,301*48);

% 40 is numOfGroups for groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar';
allSesActivityArr1 = cell(40,1);
allSesActivityArr2 = cell(40,1);

%% If data is already prepared for analysis, load the existing collected data from hard drive based on the previos section; 
% it runs way faster if you load and analyze the two arrays separately
% if data is already saved sapartely for arrays, set the path bellow where to load the data from  

if dataIsAlreadyPrep
    switch [Monkey ', ' selectedArray]
        case 'MacDuff, Anterior'
            load '';
        case 'MacDuff, Posterior'
            load '';
        case 'Mojo, Anterior'
            load '';
        case 'Mojo, Posterior'
            load '';
    end
end
%% load each channel-day file; zscore LFP within freq.; group (based on groupingCriteria);  and collect LFP from all sessions 
% if you have alreday collected and saved the data, do not run this section.
if ~dataIsAlreadyPrep

for ses = 1:length(all_clean_folder_names)
    tempActivity = cell(40,1);  % 36 is numOfGroups for groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar';
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
                % assign electrode position on the array for current channel
                eval(['plot_indx(2,ch)   = find(' 'Array' num2str(arrayIdx) ' == ch);']);
                
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
                % assign electrode position on the array for current channel
                eval(['plot_indx(2,ch)   = find(' 'Array' num2str(arrayIdx) ' == str2num(freq.label{1}(k+1:end)));']);
                clear k
       end
        
        % assign array index to this channel
        plot_indx(1,ch) = arrayIdx;
               
        
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
                
% % % %                 % ---- Event timing outlier removal: remove trials for which event times or reaction times were outliers
% % % %                 evtTimeAndRTOutTrials = cellfun(@(x,y,z) (x==1)||(y==1), {event.OutCueOnOffZ2}', {event.OutAqToRwdZ2}');
% % % %                 % remove outlier trials from event
% % % %                 event(evtTimeAndRTOutTrials) = [];
% % % %                 % remove outlier trials from LFP power of the current channel
% % % %                 Z_freq_powspctrm(:,:,evtTimeAndRTOutTrials) = [];
% % % %                 display([num2str(length(find(evtTimeAndRTOutTrials))) ' trials were removed due to bad event timing'])
% % % %                 clear evtTimeAndRTOutTrials
                
                % in order to match the task conditions of the Uncertainty task with the
                % Risk task, only the following task conditions are kept:
                % C0, C50, C100, R50
                aa = struct2table(event);
                bb = aa.type;
                c1 = cellfun(@(x) ~isempty(strfind(x,'C0')),bb,'UniformOutput',0);
                c2 = cellfun(@(x) ~isempty(strfind(x,'C50')),bb,'UniformOutput',0);
                c3 = cellfun(@(x) ~isempty(strfind(x,'C100')),bb,'UniformOutput',0);
                c4 = cellfun(@(x) ~isempty(strfind(x,'R50')),bb,'UniformOutput',0);
                ind = sort([find(cell2mat(c1));find(cell2mat(c2));find(cell2mat(c3));find(cell2mat(c4))]);
                new_event = table2struct(aa(ind,:));
                clear aa bb c1 c2 c3 c4 ind
                
                % group remaining trials based on the groupingCriteria
                [output] = GroupTrials(new_event, groupingCriteria );% 'event' only includes successful trials, so it is not necessary to remove error trials
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
eval(['allSesActivityArr' num2str(arrayIdx) ' = cell(40,1);'])

for ses = 1:length(fNames)
    % load the session data, this will load a variable named temp which is a cell array of size 36x1
    load([intermediate_save_folder   fNames(ses).name])
    for gr=1:length(temp)
        if size(allSesActivityArr2{gr},1)<40
            allSesActivityArr2{gr} = zeros(47,301,2);
        else
            
            if size(temp{gr},1)<40
                temp{gr} = zeros(47,301,2);
            end

        eval(['allSesActivityArr' num2str(arrayIdx) '{gr} = cat(3, allSesActivityArr' num2str(arrayIdx) '{gr}, temp{gr});'])
        end
    end
    clear gr temp
end
clear ses

% store  allSesActivityArr1 or allSesActivityArr2 on the hard drive
save([save_dir Monkey '_allSesActivityArr' num2str(arrayIdx)], ['allSesActivityArr' num2str(arrayIdx)]) 

end  % end of 'if ~dataIsAlreadyPrep ...'

%% Apply GLM; it runs way faster if you load and analyze the two arrays separately
% load stored power spectrum file of one channel (any channel from any session) to regenerate the groups and their parameter values
tmp_freq_folder_list = dir(freq_data_dir);
freq_file_names = tmp_freq_folder_list(3).name;
clear tmp_freq_folder_list
load ([freq_data_dir  freq_file_names])

% get the parameter values for each group
event = recursive_read(freq,'event');
aa = struct2table(event);
bb = aa.type;
c1 = cellfun(@(x) ~isempty(strfind(x,'C0')),bb,'UniformOutput',0);
c2 = cellfun(@(x) ~isempty(strfind(x,'C50')),bb,'UniformOutput',0);
c3 = cellfun(@(x) ~isempty(strfind(x,'C100')),bb,'UniformOutput',0);
c4 = cellfun(@(x) ~isempty(strfind(x,'R50')),bb,'UniformOutput',0);
ind = sort([find(cell2mat(c1));find(cell2mat(c2));find(cell2mat(c3));find(cell2mat(c4))]);
new_event = table2struct(aa(ind,:));
clear aa bb c1 c2 c3 c4 ind
[output] = GroupTrials(new_event, groupingCriteria);
% do not clear freq. it is needed bellow
clear event


        % zscore GLM predictors, get them from output and add a new field .ZValue to the output
        aa = cell2mat({output.Value}');
        [C1,~,ic1] = unique(aa(:,1));
        z1= zscore(C1);
        [C2,~,ic2] = unique(aa(:,2));
        z2= zscore(C2);
        [C3,~,ic3] = unique(aa(:,3));
        z3= zscore(C3);
        [C4,~,ic4] = unique(aa(:,4));
        z4= zscore(C4);        
        for i=1:length(ic1)
          output(i).ZValue =[z1(ic1(i))  z2(ic2(i)) z3(ic3(i)) z4(ic4(i))]; 
        end
        clear aa C1 C2 C3 C4 ic1 ic2 ic3 ic4 z1 z2 z3 z4 i

% % % % zscore regressors
% % % cuePosZ   =  [-1           1
% % %               -0.7071   0.7071];
% % % trgetPosZ =  [-1           1
% % %               -0.7071   0.7071]; 
% % % evZ       =  [3   6  9
% % %               -1  0  1];
% % % varZ      =  [0        1      4
% % %               -0.80  -0.32  1.12];
% % % for gr = 1:length(output)
% % %     idx1 = find( output(gr).Value(1) == cuePosZ(1,:) );
% % %     idx2 = find( output(gr).Value(2) == trgetPosZ(1,:) );
% % %     idx3 = find( output(gr).Value(3) == evZ(1,:) );
% % %     idx4 = find( output(gr).Value(4) == varZ(1,:) );
% % %     output(gr).Value = [cuePosZ(2,idx1), trgetPosZ(2,idx2), evZ(2,idx3), varZ(2,idx4)];
% % %     clear idx1 idx2 idx3 idx4
% % % end
% % % clear gr
    
% ---------------------------------------------------------------------------------------------------------------------------
fq = 0;
b=cell(1);
stats=cell(1); 
eval(['numFreq = size(allSesActivityArr' num2str(arrayIdx) '{1},1);'])
eval(['numTime = size(allSesActivityArr' num2str(arrayIdx) '{1},2);'])
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
            eval(['glmData = [glmData; squeeze(allSesActivityArr' num2str(arrayIdx) '{gr}(fq,t,:))];'])  % some trials have NaN value which means those electrodes must be removed because of artifact
            
            eval(['glmPred = [glmPred; repmat(output(gr).ZValue, length(squeeze(allSesActivityArr' num2str(arrayIdx) '{gr}(fq,t,:))), 1)];'])
        end
        clear gr    

%         [b{fq,t},~,stats{fq,t}] = glmfit(glmPred, glmData-min(glmData)+eps, 'gamma');
        [b{fq,t},~,stats{fq,t}] = glmfit(glmPred, glmData, 'normal');
        [fq t]
    end
end


%% Plot the GLM coefficient for each time and frequency; use False Discovery Rate for significance level

paraNames = {'CuePos', 'TargetPosition', 'PreGood', 'CurrEV'};
fg1= figure('Name',Array_pos{1}); fg2= figure('Name',Array_pos{2});
colorRange = [-0.12, 0.08];

for para = 2:5 % for each parameter except for the constant value!
    
    temp = cellfun(@(x) x(para), b, 'UniformOutput', 0);
    % get the p value for this parameter
    pVal = cellfun(@(x) x.p(para), stats, 'UniformOutput', 0); 
    pValPara = cell2mat(pVal);
    pThr = cellfun(@(x) fdr(x), mat2cell(pValPara, ones(size(pValPara,1),1), size(pValPara,2)), 'UniformOutput', 0); 
% % % %     pThr = cellfun(@(x) fdr(x.p), stats, 'UniformOutput', 0); % find the significant pThreshold defined by False Discovery Rate method
% % % %     signifP = cell2mat(pVal) <= cell2mat(pThr);
    pThr = repmat(cell2mat(pThr), 1 ,  size(pValPara,2) );
    signifP = cell2mat(pVal) <= pThr;
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
