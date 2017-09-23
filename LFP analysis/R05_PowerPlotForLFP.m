
% this script load the existing preprocessed LFP for electrod-days; zscores
% each electrode-day data within each frequency; groups trials based on groupingCriteria
% plots the LFP power for each group

clear
close all
clc

%%

Monkey            = 'Mojo';
groupingCriteria  = 'expected_reward';

% color_axis = [-0.6,1.8];% MacDuff
color_axis = [-0.5,0.7];% Mojo
num_ch = 96;

%%

if strcmp(Monkey,'MacDuff')
    % we have the channel numbers (not the electrode numbers)
    [ Array1,Array2 ] = electrodepinoutMacDuff('chan');% Array1 is anterior and Array2 is posterior
    plt_row = 10;
    plt_col = 5;
    main_folder   = 'Z:\';
    data_folder   = 'data\Risk\';
    bhv_data_dir  = 'Behavior\';
    freq_data_dir = 'Frequency Data\';
    save_plot_dir = 'D:\shared\New_TimeFrequency plot\';
elseif strcmp(Monkey,'Mojo')
    % we have the channel numbers (not the electrode numbers)
    [ Array1,Array2 ] = electrodepinoutMojo('chan');% Array1 is anterior and Array2 is posterior
    plt_row = 7;
    plt_col = 7;
    main_folder   = 'F:\';
    data_folder   = '';
    bhv_data_dir  = 'Behavior\';
    freq_data_dir = 'Frequency Data\';
    save_plot_dir = 'D:\shared\New_TimeFrequency plot\';
end

%%

Array_pos = {'Anterior','Posterior'};% Array1 is anterior and Array2 is posterior

tmp_freq_folder_list = dir([main_folder, data_folder, freq_data_dir]);
all_freq_file_names = {tmp_freq_folder_list.name};
all_freq_file_names(1:2) = [];

tmp_bhv_folder_list = dir(['Z:\', 'data\Risk\', bhv_data_dir]);
all_bhv_folder_names = {tmp_bhv_folder_list.name};
all_bhv_folder_names(1:2) = [];

% file_names = all_file_names(~cellfun(@isempty, strfind(all_file_names, Alignment)));
all_bhv_folder_names = all_bhv_folder_names(~cellfun(@isempty, strfind(all_bhv_folder_names, Monkey)));          
ses_num = 0;
Ses_sum1 = zeros(47*3,1501*48);
Ses_sum2 = zeros(47*3,1501*48);
for ses = 1:length(all_bhv_folder_names) 
    file_names = all_freq_file_names(~cellfun(@isempty, strfind(all_freq_file_names, all_bhv_folder_names{ses}(5:end-4))));
    if isempty(file_names)
        continue
    end
    ses_num = ses_num+1;
%     elec_day_array1 = {};
%     elec_day_array2 = {};    
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
        
        plot_indx(1,ch) = str2double(arrayIdx);
        eval(['plot_indx(2,ch)   = find(' 'Array' arrayIdx ' == true_ch_num(ch));']);

        
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
                    eval(['elec_day_array' arrayIdx '{gr,ch} = squeeze(mean(Z_freq_powspctrm(:,:,idx),3));'])
                end
                clear gr output event
        end  
        clear Z_freq_powspctrm Z_M M M0
        
    end
    
    if strcmp(groupingCriteria,'expected_reward')
        out = [3,6,9];
    elseif strcmp(groupingCriteria,'Variance')
        out = [0,1,4];
    end
    
    scr = get(0,'screensize');
    
    for gr = 1:3
        for ch = 1:num_ch
            ch
            H = figure(plot_indx(1,ch));
            H.Name = ['      ' groupingCriteria '= ' num2str(out(gr)) '               Array position: ' Array_pos{plot_indx(1,ch)}];
            H.Position = scr;
            
%             title([num2str(true_ch_num(ch))])
            subplot(plt_row,plt_col,plot_indx(2,ch))
            eval(['H = mesh(elec_day_array' num2str(plot_indx(1,ch)) '{gr,ch});']);
            colormap('jet'),view([0 90]),grid on,axis('tight')
            H.XData = freq.time;
            line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
            ax = gca;ax.XTickLabel(find(strcmp(ax.XTickLabel,'0'))) = {'Cue On'};
            caxis(color_axis)
            axis off
            clear H ax
            
%             eval(['MAT = elec_day_array' num2str(plot_indx(1,ch)) '{gr,ch};']);
%             Extr(ch,:) = [min(MAT(:)),max(MAT(:))];
              
        end
        figure(1)
        saveas(gcf,[save_plot_dir,all_bhv_folder_names{ses}(5:end-4),' - ',Array_pos{1},' - ',groupingCriteria,' ',num2str(out(gr)),'.jpg'])
        figure(2)
        saveas(gcf,[save_plot_dir,all_bhv_folder_names{ses}(5:end-4),' - ',Array_pos{2},' - ',groupingCriteria,' ',num2str(out(gr)),'.jpg'])
        close all

        Elec_sum_array1{gr,ses} = squeeze(mean(reshape(cell2mat(elec_day_array1(gr,:)),numel(freq.freq),length(freq.time),[]),3));% mean across electrodes
        figure
        H = mesh(Elec_sum_array1{gr,ses});
        colormap('jet'),view([0 90]),grid on,axis('tight')
        H.XData = freq.time;line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        ax = gca;ax.XTickLabel(find(strcmp(ax.XTickLabel,'0'))) = {'Cue On'};
        xlabel('Time(sec)'),ylabel('Frequency(Hz)')
        saveas(gcf,[save_plot_dir,all_bhv_folder_names{ses}(5:end-4),' - ','Array ',num2str(1),' - ',groupingCriteria,' ',num2str(out(gr)),' Electrode Averaged','.jpg'])
        clear H ax
        close all
        
        Elec_sum_array2{gr,ses} = squeeze(mean(reshape(cell2mat(elec_day_array2(gr,:)),numel(freq.freq),length(freq.time),[]),3));% mean across electrodes
        figure
        H = mesh(Elec_sum_array2{gr,ses});
        colormap('jet'),view([0 90]),grid on,axis('tight')
        H.XData = freq.time;line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        ax = gca;ax.XTickLabel(find(strcmp(ax.XTickLabel,'0'))) = {'Cue On'};
        xlabel('Time(sec)'),ylabel('Frequency(Hz)')
        saveas(gcf,[save_plot_dir,all_bhv_folder_names{ses}(5:end-4),' - ','Array ',num2str(2),' - ',groupingCriteria,' ',num2str(out(gr)),' Electrode Averaged','.jpg'])
        clear H ax
        close all
        
%         axis on
%         subplot(10,5,1),xlabel('Time(sec)'),ylabel('Frequency(Hz)'),colorbar
        
    end
    
    
    
    %% mean across sessions
    Ses_sum1 = Ses_sum1 + cell2mat(elec_day_array1);
    Ses_sum2 = Ses_sum2 + cell2mat(elec_day_array1);
    
%     if strcmp(Grouping,'EV')
%     out = [3,6,9];
% elseif strcmp(Grouping,'Var')
%     out = [0,1,4];
% end
% 
% for gr = 1:3
%         if strcmp(Grouping,'No')
%             H = mesh(norm_Ant{1});
%             colormap('jet'),colorbar,view([0 90]),grid on,axis('tight')
%             H.XData = freq.time;xlabel('Time(sec)'),ylabel('Frequency(Hz)')
%             line([0,0],[H.YData(1),H.YData(end)],'LineStyle','--','Color',[0,0,0],'LineWidth',3)
%             ax = gca;ax.XTickLabel(2) = {'Cue On'};
%             saveas(gcf,[save_plot_dir,'All Sessions',' - ','Anterior - ','No Grouping',' ',Monkey,'.jpg'])
%             clear H ax
%             H = mesh(norm_Pos{1});
%             colormap('jet'),colorbar,view([0 90]),grid on,axis('tight')
%             H.XData = freq.time;xlabel('Time(sec)'),ylabel('Frequency(Hz)')
%             line([0,0],[H.YData(1),H.YData(end)],'LineStyle','--','Color',[0,0,0],'LineWidth',3)
%             ax = gca;ax.XTickLabel(2) = {'Cue On'};
%             saveas(gcf,[save_plot_dir,'All Sessions',' - ','Posterior - ','No Grouping',' ',Monkey,'.jpg'])
%             clear H ax
%         else
%             H = mesh(elec_day_array1{gr});
%             colormap('jet'),colorbar,view([0 90]),grid on,axis('tight')
%             H.XData = freq.time;xlabel('Time(sec)'),ylabel('Frequency(Hz)')
%             line([0,0],[H.YData(1),H.YData(end)],'LineStyle','--','Color',[0,0,0],'LineWidth',3)
%             ax = gca;ax.XTickLabel(2) = {'Cue On'};
%             saveas(gcf,[save_plot_dir,'All Sessions',' - ','Anterior - ',Grouping,'= ',num2str(out(gr)),' ',Monkey,'.jpg'])
%             clear H ax
%             H = mesh(elec_day_array2{gr});
%             colormap('jet'),colorbar,view([0 90]),grid on,axis('tight')
%             H.XData = freq.time;xlabel('Time(sec)'),ylabel('Frequency(Hz)')
%             line([0,0],[H.YData(1),H.YData(end)],'LineStyle','--','Color',[0,0,0],'LineWidth',3)
%             ax = gca;ax.XTickLabel(2) = {'Cue On'};
%             saveas(gcf,[save_plot_dir,'All Sessions',' - ','Posterior - ',Grouping,'= ',num2str(out(gr)),' ',Monkey,'.jpg'])
%             clear H ax
%         end
% end



end



close all
figure
for gr = 1:3
    Ses_avg_gr = Ses_sum1((gr-1)*47+1:gr*47,:)/ses_num;
    for ch = 1:48
        Ses_avg_elec = Ses_avg_gr(:,(ch-1)*1501+1:ch*1501);
        
        H = mesh(Ses_avg_elec);
        colormap('jet'),view([0 90]),grid on,axis('tight')
        H.XData = freq.time;line([0,0],[H.YData(1),H.YData(end)],[2,2],'LineStyle',':','Color',[0,0,0],'LineWidth',2)
        ax = gca;ax.XTickLabel(find(strcmp(ax.XTickLabel,'0'))) = {'Cue On'};
        xlabel('Time(sec)'),ylabel('Frequency(Hz)')
        saveas(gcf,[save_plot_dir,all_bhv_folder_names{ses}(5:end-4),' - ','Array ',num2str(2),' - ',groupingCriteria,' ',num2str(out(gr)),' Session Averaged','.jpg'])
        clear H ax
        close all
    end
end





