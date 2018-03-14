




% last update: 04.03.2018

%%

clear all;
close all;
clc;
dbstop if error 

spike_dir           = 'X:\data\Cut and Aligned Spk\';
save_dir            = 'X:\Report 13\Spike Sorting\';


%% 

% % % % ftext = fileread(['X:\data\Original Data_Extracted\MacDuff_2015_08_05\ADJ2015_08_05MacDuff1Risk002onoff']);
% % % % 
% % % % tmp_onoff       = regexp(ftext, '{(?<on>.*)}.*{(?<off>.*)}', 'names');
% % % % onoff           = {round(cellfun(@str2num, regexp(tmp_onoff.on, '[0-9]+', 'match'))), ...
% % % %                    round(cellfun(@str2num, regexp(tmp_onoff.off, '[0-9]+', 'match')))};
% % % % 
% % % % onoff = cell2mat(onoff')';
% % % % 
% % % % clear tmp_onoff ftext
% % % % 
% % % % Expression1 = Expression1(find(Expression1(:,1)), :);
% % % % grpIDs = unique(Expression1(:,2));
% % % % 
% % % % aa = cellfun( @(x) find((onoff(:,1)<=x) .* (x<onoff(:,2))) , mat2cell(Expression1(:,3),ones(size(Expression1,1),1),1) , 'UniformOutput' , 0);


%% 
if 1
    
sessionList = dir([spike_dir]); % Bahareh: this is the list of all field data files
sessionList = sessionList(3:end);

%% 

allSesUnitStability = {[]};
allSesNumSuccessTrl = {};
for ses = 1:numel(sessionList) % for each session
    ses
    % load aligned and trial-sorted spikes for this session; this will load spikeTrial variable; 
    load([spike_dir sessionList(ses).name ] );
    
    % get the strat and stop time of this recording session
    allSesNumSuccessTrl{ses} = size(spikeTrials.cfg.event,1); 
    
    temp = cellfun(@(x) [min(x) max(x)], spikeTrials.trial, 'UniformOutput', 0)'
    temp( cellfun(@(x) isempty(x), temp) ) = [];
    allSesUnitStability{ses} = cell2mat( temp );
    clear temp
    
    fg = figure('Name', sessionList(ses).name), hold on, box on
    line([ repmat(1:size(allSesUnitStability{ses},1) , 2,1) ] , allSesUnitStability{ses}')
    set(gca, 'YLim', [1 allSesNumSuccessTrl{ses}+30])
    
    title(regexprep(sessionList(ses).name, '_', '-'))
    xlabel('unit number')
    ylabel('number of successful trials')
    
% % % %     saveas(fg, [save_dir sessionList(ses).name(1:end-4)], 'fig')
    saveas(fg, [save_dir sessionList(ses).name(1:end-4)], 'bmp')
    
    close(fg)
    clear spikeTrials fg
end
clear ses 


end
