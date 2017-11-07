


clear
close all
clc

%%

main_folder    = 'Z:\';
data_folder    = 'data\';
clean_data_dir = 'Clean Data\';

Monkey         = 'Mojo';
Alignment      = 'cue';

cfg             = [];
cfg.method      = 'wavelet';
cfg.output      = 'pow';
cfg.foi         = 4:50;
cfg.polyremoval = 1; % Linear trend
cfg.width       = realpow(exp(1), sort(linspace(log(5), log(14), 50), 'ascend'));
cfg.gwidth      = 3;
cfg.padtype     = 'zero';
cfg.pad         = 'nextpow2';
cfg.verbose     = 0;
% following is our time of interest for spectrum analysis
cfg.toi         = - 0.3: 0.001: 1.2;  
cfg.keeptrials  = 'yes';
    
%%
tmp_folder_list = dir([main_folder, data_folder, clean_data_dir]);
all_file_names = {tmp_folder_list.name};
file_names = all_file_names(~cellfun(@isempty, strfind(all_file_names, Alignment)));
file_names = file_names(~cellfun(@isempty, strfind(file_names, Monkey)));          

for i = 1: length(file_names) % for each session
    try
        % load the preprocessed LFP data together with event info
        load([main_folder, data_folder, clean_data_dir, file_names{i}]);
        
        for ch = 1: length(data.label) % for each channel
            cfg.channel     = data.label(ch);
            % get the power spectrum of the current channel for this session
            [freq] = ft_freqanalysis(cfg, data);

            event = recursive_read(freq, 'event');
            
            LFP_outliers = cellfun(@(x) unique(x'), {event.bad_lfp_channels},'UniformOutput',0)';
            
            % get the bad trial indices for the current channel
            thisChBadTrls   = find(cellfun(@(x) ~isempty(find(x == ch, 1)), LFP_outliers));
            
            % replace the powerspectrum of the bad with NaN 
            freq.powspctrm(thisChBadTrls, 1, :, :) = nan( length(thisChBadTrls), size(freq.powspctrm,3), size(freq.powspctrm,4) );
            
            % store the power spectrum of the current channel on the hard drive
            save(['E:\Electrode-day power spectrum\' Monkey '\',file_names{i}(1:end-4),' ch',num2str(ch)],'freq')
        
            % free some memory
            clear freq event thisChBadTrls
        end
    catch err
        warning([err.message ' in file: ' file_names{i}]);
        continue;
    end
    clear data
end

