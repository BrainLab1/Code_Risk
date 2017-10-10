
% this script load the existing preprocessed LFP for electrod-days; zscores
% each electrode-day data within each frequency; groups trials based on 'CuePos & TargetPos & CurrEV & CurrVar'
% the session and thefor each trial take the mean LFP power across electrodes within the session so that for each trial,
% we will have two activity map (one for each array); data from the two arrays are analyzed separately.
% collects the data from all sessions and applies GLM on it. the plots the
% GLM coefficients in the time-freq. space

clear
close all
clc

%%

Monkey            = 'Mojo';
groupingCriteria  = 'CuePos & TargetPos & CurrEV & CurrVar';

%%

if strcmp(Monkey,'MacDuff')
    % we have the channel numbers (not the electrode numbers)
    [ Array1,Array2 ] = electrodepinoutMacDuff('chan');% Array1 is anterior and Array2 is posterior
    main_folder   = 'Z:\';
    data_folder   = 'data\Risk\';
    bhv_data_dir  = 'Behavior\';
    freq_data_dir = 'Frequency Data\';
    save_dir = 'X:\Report 5\';
elseif strcmp(Monkey,'Mojo')
    % we have the channel numbers (not the electrode numbers)
    [ Array1,Array2 ] = electrodepinoutMojo('chan');% Array1 is anterior and Array2 is posterior
    main_folder   = 'F:\';
    data_folder   = '';
    bhv_data_dir  = 'Behavior\';
    freq_data_dir = 'Frequency Data\';
    save_dir = 'X:\Report 5\';
end

%%

Array_pos = {'Anterior','Posterior'};% Array1 is anterior and Array2 is posterior

% tmp_freq_folder_list = dir([main_folder, data_folder, freq_data_dir]);
% all_freq_file_names = {tmp_freq_folder_list.name};
% all_freq_file_names(1:2) = [];
% 
% tmp_bhv_folder_list = dir(['Z:\', 'data\Risk\', bhv_data_dir]);
% all_bhv_folder_names = {tmp_bhv_folder_list.name};
% all_bhv_folder_names(1:2) = [];
% 
% % file_names = all_file_names(~cellfun(@isempty, strfind(all_file_names, Alignment)));
% all_bhv_folder_names = all_bhv_folder_names(~cellfun(@isempty, strfind(all_bhv_folder_names, Monkey)));          
ses_num = 0;
Ses_sum1 = zeros(47*3,1501*48);
Ses_sum2 = zeros(47*3,1501*48);

allSesActivityArr1 = cell(36,1);
allSesActivityArr2 = cell(36,1);

%% load each channel-day file; zscore LFP within freq.; group (based on groupingCriteria);  and collect LFP from all sessions 
% if you have alreday collected and saved the data, do not run this section
% and just load the collected data from hard drive
if 0 

for ses = 1:length(all_bhv_folder_names) 
    ses
    tempActivity = cell(36,1);
    file_names = all_freq_file_names(~cellfun(@isempty, strfind(all_freq_file_names, all_bhv_folder_names{ses}(5:end-4))));
    if isempty(file_names)
        continue
    end
    ses_num = ses_num+1;
    elec_day_array1 = cell(3,96);
    elec_day_array2 = cell(3,96);
    ch = 0;
    while ch <length(file_names)
        ch = ch+1;
               
        % load variable freq; size(squeeze(freq.powspctrm)) is [numTrials, numFreqBands, numTimeStamps] 
        load([main_folder, data_folder, freq_data_dir, file_names{ch}]);
        
        if strcmp(Monkey,'Mojo')
            true_ch_num(ch) = str2double(freq.label{1}(5:end));
            if ismember(true_ch_num(ch),Array1)
                arrayIdx        = num2str(1);  % elec1:Anterior; elec2:Posterior
            else
                arrayIdx        = num2str(2);
            end
            
        else
            arrayIdx        = num2str(freq.label{1}(5));  % elec1:Anterior; elec2:Posterior
            true_ch_num(ch) = str2double(file_names{ch}(strfind(file_names{ch},'ch')+2:strfind(file_names{ch},'.mat')-1));
        end
        
        % change order of freq.powspctrm dimensions and change it to [numFreqBands, numTimeStamps, numTrials]
        M0 = permute(squeeze(freq.powspctrm(:,1,:,:)), [2,3,1]);
        % turn the M0 matrix to a 2D matrix of size [numFreqBands (numTrials*numTimaStamps)]
        M = reshape(M0, size(M0,1), []);
        % zscore ezch frequency band separately 
        Z_M = zscore(M')';
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
                    % take the mean activity map for this electrode in the current session
                    eval(['elec_day_array' arrayIdx '{gr,ch} = squeeze(mean(Z_freq_powspctrm(:,:,idx),3));'])
                end
                clear gr output event
        end  
        clear Z_freq_powspctrm Z_M M M0
        
    end
    
    % in this session, for each trial, take the mean activity across electrodes 
    arr1Idx = ismember(true_ch_num, Array1)';
    temp1 = cellfun(@(x) squeeze( mean(x(find(arr1Idx),:,:,:), 1) ),tempActivity, 'UniformOutput', false);
    arr2Idx = ismember(true_ch_num, Array2)';
    temp2 = cellfun(@(x) squeeze( mean(x(find(arr2Idx),:,:,:), 1) ),tempActivity, 'UniformOutput', false);
    
    for gr = 1:length(allSesActivityArr1)
        allSesActivityArr1{gr} = cat(3, allSesActivityArr1{gr}, temp1{gr});
        allSesActivityArr2{gr} = cat(3, allSesActivityArr2{gr}, temp2{gr});
    end
    clear gr temp tempActivity true_ch_num temp1 temp2 arr1Idx arr2Idx
    
    
end
clear ses    
    
save([save_dir Monkey '_allSesActivityArr1\'],allSesActivityArr1)    
save([save_dir Monkey '_allSesActivityArr2\'],allSesActivityArr2)

end

%% Load the existing collected data from hard drive based on the previos section; it runs way faster if you load and analyze the two arrays separately
if 1
% % % % load 'X:\Report 5\MacDuff_allSesActivityArr1';
% % % % load 'X:\Report 5\MacDuff_allSesActivityArr2';

% load 'D:\Shared\Report 5\Mojo_allSesActivityArr1';
load 'D:\Shared\Report 5\Mojo_allSesActivityArr2';
load 'Z:\data\Risk\Frequency Data\Data_MacDuff_2015_08_05 (cue) ch1';

end
%% Apply GLM; it runs way faster if you load and analyze the two arrays separately
% get the parameter values for each group
event = recursive_read(freq,'event');
[output] = GroupTrials(event, groupingCriteria);
clear event

fq = 0;
b1=cell(1); b2=cell(1);
stats1=cell(1); stats2=cell(1);
while fq < size(allSesActivityArr2{1},1) % for each frequency
    fq = fq+1;
    t = 0;
    while t < size(allSesActivityArr2{1},2) % for each time sample
        t = t+1;
        % collect data and predictor values for GLM, separately for arrays
        gr = 0;
        glmData1 = [];
%         glmData2 = [];
        glmPred1 = [];
%         glmPred2 = [];
        while gr<length(output)
            gr = gr+1;
            glmData1 = [glmData1; squeeze(allSesActivityArr2{gr}(fq,t,:))];
%             glmData2 = [glmData2; squeeze(allSesActivityArr2{gr}(fq,t,:))];
            
            glmPred1 = [glmPred1; repmat(output(gr).Value, length(squeeze(allSesActivityArr2{gr}(fq,t,:))), 1)];
%             glmPred2 = [glmPred2; repmat(output(gr).Value, length(squeeze(allSesActivityArr2{gr}(fq,t,:))), 1)];
        end
        clear gr    

%         [b1{fq,t},~,stats1{fq,t}] = glmfit(glmPred1, glmData1-min(glmData1)+eps, 'gamma');
%         [b2{fq,t},~,stats2{fq,t}] = glmfit(glmPred2, glmData2-min(glmData2)+eps, 'gamma');
        [b1{fq,t},~,stats1{fq,t}] = glmfit(glmPred1, glmData1, 'normal');
%         [b2{fq,t},~,stats2{fq,t}] = glmfit(glmPred2, glmData2, 'normal');
        [fq t]
    end
end


%% Plot the GLM coefficient for each time and frequency; use False Discovery Rate for significance level

paraNames = {'CuePos', 'TargetPosition', 'CurrEV', 'CurrVar'};
fg1= figure('Name',Array_pos{1}); fg2= figure('Name',Array_pos{2});
colorRange = [-0.12, 0.08];

for para = 2:5
    
    temp1 = cellfun(@(x) x(para), b1, 'UniformOutput', 0);
    % get the p value for this parameter
    pVal1 = cellfun(@(x) x.p(para), stats1, 'UniformOutput', 0); 
    pThr1 = cellfun(@(x) fdr(x.p), stats1, 'UniformOutput', 0); % find the significant pThreshold defined by False Discovery Rate method
    signifP1 = cell2mat(pVal1) <= cell2mat(pThr1);
    figure(fg1)
    subplot(2,2,para-1), hold on, box on
    title(paraNames{para-1})
    H = mesh(cell2mat(temp1));
    colormap('jet'),view([0 90]),grid on,axis('tight')
    caxis([min(cell2mat(temp1(:)))-0.01 max(cell2mat(temp1(:)))+0.01]);
    H.XData = freq.time;
    set(gca,'XTick', [0 0.4 0.8], 'XTickLabel', {'CueOn','CueOff','TargetOn'})
    line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    line([0.4,0.4],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    line([0.8,0.8],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    ylabel('Frequency')
%     [min(cell2mat(temp1(:))) max(cell2mat(temp1(:)))]
    colorbar
%     contour(H.XData, H.YData, signifP1, 1)
    contour(H.XData, H.YData, signifP1, 1, 'LineColor', 'k', 'ContourZLevel', 1.5);
    
%     temp2 = cellfun(@(x) x(para), b2, 'UniformOutput', 0);
%     % get the p value for this parameter
%     pVal2 = cellfun(@(x) x.p(para), stats2, 'UniformOutput', 0); 
%     pThr2 = cellfun(@(x) fdr(x.p), stats2, 'UniformOutput', 0); % find the significant pThreshold defined by False Discovery Rate method
%     signifP2 = cell2mat(pVal2) <= cell2mat(pThr2);
%     figure(fg2)
%     subplot(2,2,para-1), hold on, box on
%     title(paraNames{para-1})
%     H = mesh(cell2mat(temp2));
%     colormap('jet'),view([0 90]),grid on,axis('tight')
%     caxis(colorRange);
%     H.XData = freq.time;
%     set(gca,'XTick', [0 0.4 0.8], 'XTickLabel', {'CueOn','CueOff','TargetOn'})
%     line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
%     line([0.4,0.4],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
%     line([0.8,0.8],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
%     ylabel('Frequency')
% %     [min(cell2mat(temp2(:))) max(cell2mat(temp2(:)))]
%     colorbar
%     contour(H.XData, H.YData, signifP2, 1)
    
    
    
    clear temp1 temp2 pVal1 pThr1 signifP1 pVal2 pThr2 signifP2
end
