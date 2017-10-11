
% last modified by Bahareh; 02.10.2017; there was an error for assigning electrode location for some sessions of Mojo when channels were named
%      elec#-#, which was fixed. session power spectrums which were plotted before this date were affected by the error
% last modified: Saeed; 16.09.2017
% 13.09.2017 Saeed

clear
close all
clc

%% set some parameters
Monkey            = 'Mojo';
groupingCriteria  = 'expected_reward';
num_ch = 96;

%%
if strcmp(Monkey,'MacDuff')
    plt_row = 10;
    plt_col = 5;
    main_folder   = 'F:\';
    data_folder   = '';
    bhv_data_dir  = 'Behavior\';
    freq_data_dir = 'Electrode-day power spectrum\MacDuff\';
    save_plot_dir = 'D:\Shared\LFP power plots\';
    intermediate_save_dir = ['D:\Shared\LFP power plots\MacDuff\' groupingCriteria '\intermediate save\'];
    color_axis = [-0.6,1.8];

elseif strcmp(Monkey,'Mojo')
    plt_row = 7;
    plt_col = 7;
    main_folder   = 'F:\';
    data_folder   = '';
    bhv_data_dir  = 'Behavior\';
    freq_data_dir = 'Electrode-day power spectrum\Mojo\';
    save_plot_dir = 'D:\Shared\LFP power plots\';
    intermediate_save_dir = ['D:\Shared\LFP power plots\Mojo\' groupingCriteria '\intermediate save\'];
    color_axis = [-0.5,0.7];
end

%%

Array_pos = {'Anterior','Posterior'};% Array1 is anterior and Array2 is posterior

tmp_freq_folder_list = dir([main_folder, data_folder, freq_data_dir]);
all_freq_file_names = {tmp_freq_folder_list.name};
all_freq_file_names(1:2) = [];

tmp_bhv_folder_list = dir(['Z:\', 'data\Risk\', bhv_data_dir]);
all_bhv_folder_names = {tmp_bhv_folder_list.name};
all_bhv_folder_names(1:2) = [];

all_bhv_folder_names = all_bhv_folder_names(~cellfun(@isempty, strfind(all_bhv_folder_names, Monkey)));          
if strcmp(groupingCriteria,'All trials')
  % here is what the numbers mean: 47 is numOfFreqBands; 1501 is numOfTimeStamps; 48 is numOfElectrodsInEachArray  
    Ses_sum1 = zeros(47,1501*48); % for array1 (anterior)
    Ses_sum2 = zeros(47,1501*48); % for array2 (posterior)
else
  % here is what the numbers mean: 47 is numOfFreqBands; 3 is numOfGroups for groupingCriteria  = 'expected_reward';
  % 1501 is numOfTimeStamps; 48 is numOfElectrodsInEachArray  
    Ses_sum1 = zeros(47*3,1501*48); % for array1 (anterior)
    Ses_sum2 = zeros(47*3,1501*48); % for array2 (posterior)
end

for ses = 1:length(all_bhv_folder_names) 
    
    % read out all the file names for channels of this session
    file_names = all_freq_file_names(~cellfun(@isempty, strfind(all_freq_file_names, all_bhv_folder_names{ses}(5:end-4))));
    if isempty(file_names)
        continue
    end
    
    % sort the file names based on the channel lable index (i.e., X in the 'chX' of the freq. file name !!)
    a1 = cellfun(@(x)    str2double(x(max(strfind(x,'ch'))+2 : min(strfind(x,'.'))-1)), file_names);
    [~,sortedIdx]=sort(a1);
    file_names = file_names(sortedIdx);
    clear a1 sortedIdx
    
    if strcmp(groupingCriteria,'All trials')
        elec_day_array1 = cell(1,96);
        elec_day_array2 = cell(1,96);
    else
      % 3 represents number of groups for groupingCriteria of 'expected_reward'. it has to change if the number of groups chahge
        elec_day_array1 = cell(3,96);
        elec_day_array2 = cell(3,96);  
    end
    
	% collect the power spectrum from all channels in this session, separated by array
    ch = 0;
    while ch <length(file_names)
        ch = ch+1;
               
        % load variable freq; size(squeeze(freq.powspctrm)) is [numTrials, numFreqBands, numTimeStamps] 
        load([main_folder, data_folder, freq_data_dir, file_names{ch}]);
        
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
                        display('There is an error! This channels belong to neither of the arrays!!!')
                        break
                    end
                end
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
                eval(['plot_indx(2,ch)   = find(' 'Array' num2str(arrayIdx) ' == str2num(freq.label{1}(k+1:end)));']);
                clear k
       end
        
        % assign array index to this channel
        plot_indx(1,ch) = arrayIdx;
               
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
                eval(['elec_day_array' num2str(arrayIdx) '{1,ch} = squeeze(nanmean(Z_freq_powspctrm(:,:,:),3));'])
            otherwise
                event = recursive_read(freq,'event');
                [output] = GroupTrials(event, groupingCriteria );% 'event' only includes successful trials, so it is not necessary to remove error trials
                for gr = 1:numel(output)
                    idx = output(gr).TrialIdx;                    
                    eval(['elec_day_array' num2str(arrayIdx) '{gr,ch} = squeeze(nanmean(Z_freq_powspctrm(:,:,idx),3));'])
                end
                clear gr output event
        end  
        clear Z_freq_powspctrm Z_M M M0 arrayIdx
        
    end
   
    
                
    scr = get(0,'screensize');
    
	% Now plot the power spectrun for current session -----------------------------------------------------------------------
    switch groupingCriteria
        case 'All trials'
            numGroups = 1;
        otherwise
            numGroups = size(elec_day_array1,1);
    end
    
    % make a variable for naming the figures
    dummy = squeeze(nanmean(reshape(cell2mat(elec_day_array1), numGroups*length(freq.freq), length(freq.time), []), 3));
    colorScaleElecAve1 = [min(dummy(:)) max(dummy(:))];
    clear dummy 
    dummy = squeeze(nanmean(reshape(cell2mat(elec_day_array2), numGroups*length(freq.freq), length(freq.time), []), 3));
    colorScaleElecAve2 = [min(dummy(:)) max(dummy(:))];
    clear dummy 
    
    switch groupingCriteria
        case 'expected_reward'
            out = {'3','6','9'};
        case 'Variance'
            out = {'0','1','4'};
        case 'All trials'
            out = ' ' ;
    end
    
    for gr = 1:numGroups
	    % plot the power spectrums, for each channel separately, map the channel location on the array 
        F1 = figure('Name', [groupingCriteria ' ' out{gr} ', Array Position:' Array_pos{1}], 'Position',  scr);
        F2 = figure('Name', [groupingCriteria ' ' out{gr} ', Array position:' Array_pos{2}], 'Position',  scr);
        for ch = 1:num_ch
            ch
            eval(['figure(F' num2str(plot_indx(1,ch)) ')'])
            subplot(plt_row,plt_col,plot_indx(2,ch))
            eval(['H' num2str(plot_indx(1,ch)) ' = mesh(elec_day_array' num2str(plot_indx(1,ch)) '{gr,ch});']);
            colormap('jet'),view([0 90]), axis('tight')
            eval(['H' num2str(plot_indx(1,ch)) '.XData = freq.time;'])
            eval(['H' num2str(plot_indx(1,ch)) '.YData = freq.freq;'])
            eval(['yData = H' num2str(plot_indx(1,ch)) '.YData']);
            line([0,0],[yData(1),yData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
            text(-0.05, 10, 2, 'CueOn', 'Rotation', 90, 'FontSize', 8, 'Color', [1 1 1])
            line([0.4,0.4],[yData(1),yData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
            text(0.35, 10, 2, 'CueOff', 'Rotation', 90, 'FontSize', 8, 'Color', [1 1 1])
            line([0.8,0.8],[yData(1),yData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
            text(0.75, 10, 2, 'TargetOn', 'Rotation', 90, 'FontSize', 8, 'Color', [1 1 1])               
            caxis(color_axis)
            set(gca, 'XTick', [], 'YTick', [4 10 20 30 40], 'FontSize', 6)
            clear H* yData
        end
		
	    % Store the figures on the hard drive
        saveas(F1, [save_plot_dir, all_bhv_folder_names{ses}(5:end-4) ' - ' Array_pos{1} ' - ' groupingCriteria ' ' out{gr} '.jpg'])
        saveas(F2, [save_plot_dir, all_bhv_folder_names{ses}(5:end-4) ' - ' Array_pos{2} ' - ' groupingCriteria ' ' out{gr} '.jpg'])
        close all
        
    	% For the current session, plot averaged power spectrum across electrodes on each array
        F = figure('Name', [groupingCriteria ' ' out{gr}']);  
        % ------ plot Array 1 -------
        subplot(1,2,1), hold on
        title(['Array Position:' Array_pos{1}])
        Elec_mean_array1{gr,ses} = squeeze(nanmean(reshape(cell2mat(elec_day_array1(gr,:)), numel(freq.freq), length(freq.time), []), 3));% mean across electrodes
        H3 = mesh(Elec_mean_array1{gr,ses});
        colormap('jet'), view([0 90]), axis('tight')
        H3.XData = freq.time;
        H3.YData = freq.freq;
        line([0,0],[H3.YData(1),H3.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(-0.05, 30, 2, 'CueOn', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])
        line([0.4,0.4],[H3.YData(1),H3.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(0.35, 30, 2, 'CueOff', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])
        line([0.8,0.8],[H3.YData(1),H3.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(0.75, 30, 2, 'TargetOn', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])               
        xlabel('Time(sec)'), ylabel('Frequency(Hz)')
        caxis(colorScaleElecAve1)
        set(gca, 'XTick', [0 0.4 0.8], 'YTick', [4 10 20 30 40], 'FontSize', 10, 'PlotBoxAspectRatio', [1 1.1 1])
        
        % ------ plot Array 2 -------
        Elec_mean_array2{gr,ses} = squeeze(nanmean(reshape(cell2mat(elec_day_array2(gr,:)),numel(freq.freq),length(freq.time),[]),3));% mean across electrodes
        subplot(1,2,2), hold on
        title(['Array Position:' Array_pos{2}])
        H4 = mesh(Elec_mean_array2{gr,ses});
        colormap('jet'), view([0 90]), axis('tight')
        H4.XData = freq.time;
        H4.YData = freq.freq;
        line([0,0],[H4.YData(1),H4.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(-0.05, 30, 2, 'CueOn', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])
        line([0.4,0.4],[H4.YData(1),H4.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(0.35, 30, 2, 'CueOff', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])
        line([0.8,0.8],[H4.YData(1),H4.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(0.75, 30, 2, 'TargetOn', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])               
        xlabel('Time(sec)'), ylabel('Frequency(Hz)')
        caxis(colorScaleElecAve2)
        set(gca, 'XTick', [0 0.4 0.8], 'YTick', [4 10 20 30 40], 'FontSize', 10, 'PlotBoxAspectRatio', [1 1.1 1])

        % Save the figure on the hard drive
        saveas(F,[save_plot_dir,all_bhv_folder_names{ses}(5:end-4) ' - '  groupingCriteria,' ' out{gr} ' -ElectrodeAveraged','.jpg'])
        clear H3 H4 F
        close all
    end
    clear colorScaleElecAve1 colorScaleElecAve2
 
    % save the power spectrum of current session (Array 1 and 2) separated by condition and electrode [condition * number of electrodes]    
    save([intermediate_save_dir all_bhv_folder_names{ses}(5:end-4)], 'elec_day_array1','elec_day_array2')

    clear elec_day_array1 elec_day_array2

end

%% load the power spectrum of each session then bind them into a new variable

% read out intermediatly saved file names (each file is one session)
fName = dir(intermediate_save_dir);

% remove '.' and '..'
fName(1:2) = []; 

for f=1:length(fName)
    load([intermediate_save_dir fName(f).name])
    if ~exist('nanElec1','var')  % shows which electrodes are nan for each session and condition
        nanElec1=cell(size(elec_day_array1));    
        nanElec2=cell(size(elec_day_array2));   
    end
    % find nan electrodes for the current session
    gr=0; el=0;
    while gr<size(nanElec1,1)
        gr = gr+1;
        el=0;
        while el<size(nanElec1,2)
            el = el+1;
            if ~isempty(elec_day_array1{gr,el})
                nanElec1{gr,el} = [nanElec1{gr,el}; isnan(sum(elec_day_array1{gr,el}(:)))];
            end
            if ~isempty(elec_day_array2{gr,el})
                nanElec2{gr,el} = [nanElec2{gr,el}; isnan(sum(elec_day_array2{gr,el}(:)))];
            end
        end
    end   
    Ses_sum1 = squeeze(nansum( cat(3, Ses_sum1, cell2mat(elec_day_array1)), 3) );
    Ses_sum2 = squeeze(nansum( cat(3, Ses_sum2, cell2mat(elec_day_array2)), 3) );
    
    clear elec_day_array1 elec_day_array2
    
end


%% Plot averaged power spectrum for individual electrodes across sessions
close all

elecMapArray1 = plot_indx(2, plot_indx(1,:)==1);
elecMapArray2 = plot_indx(2, plot_indx(1,:)==2);
   
% remove other array electrodes!!
nanElec1(:, plot_indx(1,:)==2) = [];
nanElec2(:, plot_indx(1,:)==1) = [];

for gr = 1:numGroups
    F1 = figure('Name', [groupingCriteria ' ' out{gr} ', Array Position:' Array_pos{1} ', Ave across sessions'], 'Position',  scr);
    F2 = figure('Name', [groupingCriteria ' ' out{gr} ', Array position:' Array_pos{2} ', Ave across sessions'], 'Position',  scr);
    
    % extract current group from the sum across sessions
    Ses_sum_gr_array1 = Ses_sum1((gr-1)*length(freq.freq)+1 : gr*length(freq.freq), :); 
    Ses_sum_gr_array2 = Ses_sum2((gr-1)*length(freq.freq)+1 : gr*length(freq.freq), :); 
   
    for ch = 1:48
        % extract channel index 'ch' on the array and average across avalable eceltrode-sessions
        Ses_avg_elec_array1 = Ses_sum_gr_array1(:,(ch-1)*length(freq.time)+1 : ch*length(freq.time))/max(1, sum(~(nanElec1{gr,ch})));
        Ses_avg_elec_array2 = Ses_sum_gr_array2(:,(ch-1)*length(freq.time)+1 : ch*length(freq.time))/max(1, sum(~(nanElec1{gr,ch})));
        
        % plot Array1
        figure(F1)
        subplot(plt_row, plt_col, elecMapArray1(ch)), hold on, box on
        H = mesh(Ses_avg_elec_array1);
        colormap('jet'),view([0 90]), axis('tight')
        H.XData = freq.time;
        H.YData = freq.freq;
        line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(-0.05, 10, 2, 'CueOn', 'Rotation', 90, 'FontSize', 8, 'Color', [1 1 1])
        line([0.4,0.4],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(0.35, 10, 2, 'CueOff', 'Rotation', 90, 'FontSize', 8, 'Color', [1 1 1])
        line([0.8,0.8],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(0.75, 10, 2, 'TargetOn', 'Rotation', 90, 'FontSize', 8, 'Color', [1 1 1])               
        caxis(color_axis)
        set(gca, 'XTick', [], 'YTick', [4 10 20 30 40], 'FontSize', 6)
        clear H
        
        % plot Array2
        figure(F2)
        subplot(plt_row, plt_col, elecMapArray2(ch)), hold on, box on
        H = mesh(Ses_avg_elec_array2);
        colormap('jet'),view([0 90]), axis('tight')
        H.XData = freq.time;
        H.YData = freq.freq;
        line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(-0.05, 10, 2, 'CueOn', 'Rotation', 90, 'FontSize', 8, 'Color', [1 1 1])
        line([0.4,0.4],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(0.35, 10, 2, 'CueOff', 'Rotation', 90, 'FontSize', 8, 'Color', [1 1 1])
        line([0.8,0.8],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        text(0.75, 10, 2, 'TargetOn', 'Rotation', 90, 'FontSize', 8, 'Color', [1 1 1])               
        caxis(color_axis)
        set(gca, 'XTick', [], 'YTick', [4 10 20 30 40], 'FontSize', 6)
    end
    saveas(F1, [save_plot_dir 'Array ' Array_pos{1} ' - ',groupingCriteria ' ' out{gr} 'Elec Map Session Averaged','.jpg'])
    saveas(F2, [save_plot_dir 'Array ' Array_pos{2} ' - ',groupingCriteria ' ' out{gr} 'Elec Map Session Averaged','.jpg'])
    saveas(F1, [save_plot_dir 'Array ' Array_pos{1} ' - ',groupingCriteria ' ' out{gr} 'Elec Map Session Averaged','.fig'])
    saveas(F2, [save_plot_dir 'Array ' Array_pos{2} ' - ',groupingCriteria ' ' out{gr} 'Elec Map Session Averaged','.fig'])
    clear F1 F2 
    close all
end

%% Plot averaged power spectrum for each array (averaged over electrode within the session) across sessions
close all
colorScaleElecAve1 = [-0.5 1];
colorScaleElecAve2 = [-0.5 1];

for gr = 1:numGroups
    F = figure('Name', [groupingCriteria ' ' out{gr} ', Ave across sessions'], 'Position',  scr);
    
    % extract current group power for all sessions
    elecAveArray1 = cell2mat(Elec_mean_array1(gr,:)); % size = [ (numFreq) x (numTime x numSessions) ]
    elecAveArray2 = cell2mat(Elec_mean_array2(gr,:)); % size = [ (numFreq) x (numTime x numSessions) ]
    
    % take the averaged power across sessions for the two arrays
    sesAveArray1 = squeeze( nanmean(reshape(elecAveArray1, size(elecAveArray1,1), length(freq.time), []), 3) );
    sesAveArray2 = squeeze( nanmean(reshape(elecAveArray2, size(elecAveArray2,1), length(freq.time), []), 3) );

    % plot Array1
    figure(F)
    subplot(1, 2, 1), hold on, box on
    title(['Array Position: ' Array_pos{1} ])
    H = mesh(sesAveArray1);
    colormap('jet'),view([0 90]), axis('tight')
    H.XData = freq.time;
    H.YData = freq.freq;
    line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    text(-0.05, 30, 2, 'CueOn', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])
    line([0.4,0.4],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    text(0.35, 30, 2, 'CueOff', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])
    line([0.8,0.8],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    text(0.75, 30, 2, 'TargetOn', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])               
    xlabel('Time (sec)'), ylabel('Frequency (Hz)')
%     caxis(colorScaleElecAve1)
    set(gca, 'XTick', [0 0.4 0.8], 'YTick', [4 10 20 30 40], 'FontSize', 10, 'PlotBoxAspectRatio', [1 1.1 1])
    clear H
        
    % plot Array2
    subplot(1, 2, 2), hold on, box on
    title(['Array Position: ' Array_pos{2} ])
    H = mesh(sesAveArray2);
    colormap('jet'),view([0 90]), axis('tight')
    H.XData = freq.time;
    H.YData = freq.freq;
    line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    text(-0.05, 30, 2, 'CueOn', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])
    line([0.4,0.4],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    text(0.35, 30, 2, 'CueOff', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])
    line([0.8,0.8],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
    text(0.75, 30, 2, 'TargetOn', 'Rotation', 90, 'FontSize', 14, 'Color', [1 1 1])               
    xlabel('Time (sec)'), ylabel('Frequency (Hz)')
%     caxis(colorScaleElecAve2)
    set(gca, 'XTick', [0 0.4 0.8], 'YTick', [4 10 20 30 40], 'FontSize', 10, 'PlotBoxAspectRatio', [1 1.1 1])
    clear H

    saveas(F, [save_plot_dir 'Array ' Array_pos{1} ' - ',groupingCriteria ' ' out{gr} 'Session Average','.jpg'])
    saveas(F, [save_plot_dir 'Array ' Array_pos{1} ' - ',groupingCriteria ' ' out{gr} 'Session Average','.fig'])
    clear F 
    close all
end



