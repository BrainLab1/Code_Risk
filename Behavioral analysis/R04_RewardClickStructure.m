



% created on 28.08.2017 by Bahareh: This script plots the reaward signal in
% trials and shows when in a trial reward clicks happen

clear all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
data_folder            = 'Z:\data\Risk\Behavior\';

colorCode = [  204 204 255  % F3
               255 153 153  % F6
               153 255 204   % F9
               102 102 255   % B3L
               255 51  51    % B6L
               51  255 102   % B9L              
               0   0   255   % B3H
               204 0   0     % B6H
               0  153 153]/255; % B9H   
           
monkeyName = 'MacDuff';  % 'Moj' or 'Mac'
groupingCriteria       = 'expecter_reward & RewardVariance & Outcome' ; % this parameter shows based on what criteria trials are grouped.

%% add path
addpath(genpath('D:\Code_localGit'))

%% Read out list of all the files related to this session
dataPath = [data_folder];
allFiles = dir(dataPath);    

% remove the '.' and '..' from the folders list and filter sessions for monkey names
idx = [];
for ses = 1:numel(allFiles) % for each session
    if (numel(allFiles(ses).name) < 3) || ~(strcmp(allFiles(ses).name(5:7), monkeyName(1:3)))  
        idx = [idx; ses];
    end
end
allFiles(idx) = [];
clear ses idx

%% Extract and group reaction times
allSesGroupedData = {};
allSesRwdSigTarAcq = cell(15,1);  % collected reaction times from successful trials of each session

for ses = 1:numel(allFiles) % for each session
    load ([dataPath allFiles(ses).name])
    
    % convert the event tructure to table
    eventTable = struct2table(new_cfg.event);
    
    % group trials for 9 task conditions
    [output] = GroupTrials(table2struct(eventTable), groupingCriteria);

    % get the trial index for successful trials per condition
    for gr = 1:length(output)
        output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
        % for successful trials of this group, collect the RT 
        grRwSgnTarAcq = zeros(length(output(gr).SuccessTrlIdx),1100); % current group reward signals
        for tr = 1:length(output(gr).SuccessTrlIdx) % for each successful trial in this group
            trIdx = output(gr).SuccessTrlIdx(tr);

            % read out the reward timing and make the reward signal time for this successful trial of this group
            onSmpl  = eventTable.RewardOnTime{trIdx} - eventTable.ActualEventTime(trIdx).TargetAcquired;
            offSmpl = eventTable.RewardOffTime{trIdx} - eventTable.ActualEventTime(trIdx).TargetAcquired;
            
            if ~(sum(onSmpl<0) || sum(offSmpl<0))  % this line must be removed when outliers were removed. this is the 178ms exception!   
            for clk = 1:length(onSmpl)
                grRwSgnTarAcq(tr,onSmpl(clk):offSmpl(clk)) = 1;
            end
            end
            clear trIdx clk onSmpl offSmpl
        end
        output(gr).RewardSignTarAcq  = grRwSgnTarAcq; 
        allSesRwdSigTarAcq(gr) = {[allSesRwdSigTarAcq{gr}; output(gr).RewardSignTarAcq]};
        clear  grRwSgnTarAcq tr 

    end
    
    allSesGroupedData{ses} = output;
    clear gr output eventTable new_cfg 
end
clear ses

%% plot
preSetText = {'Lose', 'FixReward', 'Win'};
for gr=1:length(allSesRwdSigTarAcq)
    fg = figure, hold on, box on
    xlabel('Time (ms)')
    ylabel('Trials')
    set(fg, 'OuterPosition', [2364  508  440  246])
    mesh(allSesRwdSigTarAcq{gr, 1})
    view(2)
    yLimit = ylim;
    set(gca,'XTick',[0 377 1000],'XTickLabel',{'TargetAcquired','377ms','1000'},...
            'YTick', size(allSesRwdSigTarAcq{gr, 1},1) )
    title(['EV=' num2str(allSesGroupedData{1}(gr).Value(1)) ', Var=' num2str(allSesGroupedData{1}(gr).Value(2)) ', ' preSetText{allSesGroupedData{1}(gr).Value(3)+2}])
    saveas(fg, ['D:\group_' num2str(gr)], 'fig')
%     set(gcf, 'Renderer', 'painters')
%     saveas(fg, ['D:\group_' num2str(gr)], 'eps')
    close(fg); clear fg
end

