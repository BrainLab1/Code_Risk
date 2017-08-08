
clear all;
close all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
main_folder          = 'Z:\'; 
data_folder          = 'data\'; % 
original_data_folder = 'Risk\Original Data_Extracted\';
save_dir             = 'X:\Report 2\Bahareh\';

% Set the options
groupingCriteria       = '18 task conditions' ; % this parameter shows based on what criteria trials are grouped.

%% Read out list of all the files related to this session
dataPath = [main_folder data_folder original_data_folder];
allFiles = dir(dataPath);                                 

%%
% load matlab_fig1
if 1
    
    
numSuccess = [];
groupValue = {};
for ses = 1:numel(allFiles) % for each session
    % select the right file
    if numel(allFiles(ses).name) > 3
        sessionFolder = allFiles(ses).name;
        bhvFileInfo = dir([dataPath sessionFolder '\*.bhv']); 
        bhvFilePath = [dataPath sessionFolder '\' bhvFileInfo.name];
        
        lfpFileInfo    = dir([dataPath sessionFolder '\*.ns2']);   
        if length(lfpFileInfo) > 1
            warning(['Multiple NS files in folder ' [dataPath sessionFolder]]);
            continue;
        end       
        lfpFilePath = [dataPath sessionFolder '\' lfpFileInfo.name];

        tmp_cfg                 = [];
        tmp_cfg.headerformat    = 'blackrock_nsx';
        tmp_cfg.data_dir        = dataPath;
        tmp_cfg.trialfun        = 'fieldtrip_trialfun_RiskBhv'; % string with function name, see below (default = 'ft_trialfun_general')
        tmp_cfg.trialdef.pre    = 1; % time in seconds; this only indicates absolut value, do not use '-'
        tmp_cfg.trialdef.post   = 1; % time in seconds
        
        % Specific config for one session
        tmp_cfg.headerfile  = lfpFilePath;
        tmp_cfg.dataset     = lfpFilePath;
        tmp_cfg.session_dir = sessionFolder;
        
        tmp_cfg.trialdef.eventtype  = 'cue';
        tmp_cfg.trialdef.interval   = [-0.1 1];

        % read ou the data the same way it is done in ReadRiskDataBhv.m
        new_cfg = ft_definetrial(tmp_cfg);
        
        % group trials for 18 task conditions
        [output] = GroupTrials(new_cfg.event, groupingCriteria);

        % convert the event tructure to table
        eventTable = struct2table(new_cfg.event);
        
        % get the trial index for successful trials per condition
        for gr = 1:length(output)
            output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
            numSuccess(ses-2,gr) = length(output(gr).SuccessTrlIdx);      % size(numSuccess) = [numSessions , numGroups]
            groupValue{ses}(gr,:) = output(gr).Value;  % size(groupValue{ses}) = [numGroups x 3] -> 3 : [ExpValue Variance Congruency]
        end
        clear gr output eventTable new_cfg tmp_cfg
    end
end
clear ses


end

%% set where each condition is ploted
condPos = [{[3 0 1]}  {[6 0 1]}  {[9 0 1]};
           {[3 1 1]}  {[6 1 1]}  {[9 1 1]};
           {[3 4 1]}  {[6 4 1]}  {[9 4 1]};
           {[3 0 0]}  {[6 0 0]}  {[9 0 0]};
           {[3 1 0]}  {[6 1 0]}  {[9 1 0]};
           {[3 4 0]}  {[6 4 0]}  {[9 4 0]}];

%% plot each condition in its predefined position
% size(numSuccess)       = [numSessions , numGroups]
% size(groupValue{ses})  = [numGroups x 3] -> 3 : [ExpValue Variance Congruency]
numRows = 6;
numColumns = 3;
binWidth = 1;
xLimit = [min(min(numSuccess)) max(max(numSuccess))] + [-1 1];
yLimit = [0 15];
xbins = xLimit(1) : binWidth : xLimit(2);
figure,

for rw = 1:numRows
    for cl = 1:numColumns
        % find the frame number in the subplots
        plotIdx = cl+((rw-1)*numColumns);  
        % find which column of the numSuccess matrix should be plotted
        [~,rowNum] = ismember( condPos{rw,cl}, groupValue{3}, 'rows');
        % get the counts per bin
        [counts,~] = hist(numSuccess(:, rowNum), xbins);
        % get the max rep for this task condition
        maxRep = max(numSuccess(:,rowNum));
        % prep the plot outline
        subplot( numRows, numColumns, plotIdx ), hold on, box on
        aa = num2str(condPos{rw,cl});
        title(['EV=' aa(1) ', Var=' aa(4) ', Congruence=' aa(7)])
        % plot the lines at 10 and max  repetition for this condition
        plot([10 10] , yLimit, '-- k')
        plot([maxRep maxRep] , yLimit, '-- r')
        % plot the histogram for this condition
        bar(xbins, counts, 'FaceColor', [153 255 204]/255)
        set(gca, 'YLim', yLimit, 'XLim', xLimit, 'XTick', 0:5:xLimit(2), 'PlotBoxAspectRatio', [2, 0.5, 1])
    end
end

subplot( numRows, numColumns, 16 )
xlabel('correct trials per condition')
subplot( numRows, numColumns, 17 )
xlabel('correct trials per condition')
subplot( numRows, numColumns, 18 )
xlabel('correct trials per condition')

subplot( numRows, numColumns, 4 )
ylabel('number of sessions')
subplot( numRows, numColumns, 10 )
ylabel('number of sessions')
subplot( numRows, numColumns, 16 )
ylabel('number of sessions')









