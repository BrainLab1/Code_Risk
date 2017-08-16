
% created on 12.08.2017 by Bahareh. 
% this scriptsreads out the reaction times session-by-session; applies the
% zscore(log(reaction time)); splits the trials by lateralite of the
% saccade and plots the distributions for two monkeys. You need the
% behavioral data to have been already extracted and stored on the hard
% drive. Use RedRiskDataBhv.m to get the behavioral data in proper format
% ---------------------------------------------------------------------------------------------------------------------------

clear all;
close all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
main_folder          = 'Z:\'; 
data_folder          = 'data\'; % 
original_data_folder = 'Risk\Original Data_Extracted\';
save_dir             = 'Z:\data\Risk\Behavior\';

% Set the options
groupingCriteria       = 'SaccadeLaterality' ; % this parameter shows based on what criteria trials are grouped.

%% Read out list of all the files related to this session
dataPath = [main_folder data_folder original_data_folder];
allFiles = dir(dataPath);    

% remove the '.' and '..' from the folders list and filter sessions for monkey names
idx = [];
for ses = 1:numel(allFiles) % for each session
    if numel(allFiles(ses).name) < 3
        idx = [idx; ses];
    end
end
allFiles(idx) = [];
clear ses idx


%%

    
allSesSuccessRT = {}; % reaction time of all successful trials from all sessions and the group identifier -> [ {[reatTime saccadeTarget]} , {monkeyName}, {zscore(reactionTime)} ] 
for ses = 1:numel(allFiles) % for each session
    sessionFolder = allFiles(ses).name;
    load ([save_dir 'Bhv_' allFiles(ses).name])
        
    % group trials for 18 task conditions
    [output] = GroupTrials(new_cfg.event, groupingCriteria);

    % convert the event tructure to table
    eventTable = struct2table(new_cfg.event);
        
    % get the trial index for successful trials per condition
    temp = [];
    for gr = 1:length(output)
        output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
        temp = [temp; [eventTable.ReactionTime(output(gr).SuccessTrlIdx) , repmat( output(gr).Value, length(output(gr).SuccessTrlIdx), 1)]];
    end
    allSesSuccessRT = [allSesSuccessRT; [{temp} {allFiles(ses).name(1:3)}]];
    clear gr output eventTable new_cfg tmp_cfg temp 
end
clear ses

%% transform reaction time for each session
for ses = 1:length(allSesSuccessRT)
    allSesSuccessRT(ses,3) = {zscore( log( allSesSuccessRT{ses,1}(:,1) ) )};
end
clear ses

%% 
zScoreMacRT = [];
sesIdx = find( strcmp(allSesSuccessRT(:,2), 'Mac') );
for ses = 1:length(sesIdx)
    temp = allSesSuccessRT{sesIdx(ses),1};
    zScoreMacRT = [ zScoreMacRT; [allSesSuccessRT{sesIdx(ses),3}, temp(:,2)] ];
    clear temp
end
clear sesIdx

zScoreMojRT = [];
sesIdx = find( strcmp(allSesSuccessRT(:,2), 'Moj') );
for ses = 1:length(sesIdx)
    temp = allSesSuccessRT{sesIdx(ses),1};
    zScoreMojRT = [ zScoreMojRT; [allSesSuccessRT{sesIdx(ses),3}, temp(:,2)] ];
    clear temp
end
clear sesIdx
%% set where each condition is ploted
condPos = [{'1Moj'}   {'1Mac'};
           {'-1Moj'}  {'-1Mac'}];
       
%% plot each condition in its predefined position
numRows = 2;
numColumns = 2;
binWidth = 0.1;
xLimit = [min([zScoreMojRT(:,1);zScoreMacRT(:,1)]) , max([zScoreMojRT(:,1);zScoreMacRT(:,1)])] + [-1 1];
yLimit = [0 350];
xbins = xLimit(1) : binWidth : xLimit(2);
plotBoxAspectRatio = [2, 1, 1];

figure,

xLimit = [min(zScoreMojRT(:,1)) , max(zScoreMojRT(:,1))] + [-1 1];
% plot Mojo right saccade
subplot( numRows, numColumns, 1 ), hold on, box on
title('Mojo right saccade')
idx = find(zScoreMojRT(:,2)>0);
[counts,~] = hist(zScoreMojRT(idx,1), xbins);
bar(xbins, counts, 'FaceColor', [153 255 204]/255)
plot([0 0] , yLimit, '-- k')
plot([2 2] , yLimit, '-- r', [-2 -2] , yLimit, '-- r')
set(gca, 'YLim', yLimit, 'XLim', xLimit, 'XTick', sort(unique( [floor(xLimit(1)):2:ceil(xLimit(2)) 2 -2] )), 'PlotBoxAspectRatio', plotBoxAspectRatio)
xlabel('zscore( Log( ReactionTime ) )')
ylabel('number of samples')
dummy = [ sum(zScoreMojRT(idx,1)<-2) sum(zScoreMojRT(idx,1)>2) ];
text(-4,250, [num2str(dummy(1)) ' trials'], 'Color', 'r')
text(2.5,250, [num2str(dummy(2)) ' trials'], 'Color', 'r')
clear  idx counts dummy

% plot Mojo left saccade
subplot( numRows, numColumns, 3 ), hold on, box on
title('Mojo left saccade')
idx = find(zScoreMojRT(:,2)<0);
[counts,~] = hist(zScoreMojRT(idx,1), xbins);
bar(xbins, counts, 'FaceColor', [153 255 204]/255)
plot([0 0] , yLimit, '-- k')
plot([2 2] , yLimit, '-- r', [-2 -2] , yLimit, '-- r')
set(gca, 'YLim', yLimit, 'XLim', xLimit, 'XTick', sort(unique( [floor(xLimit(1)):2:ceil(xLimit(2)) 2 -2] )), 'PlotBoxAspectRatio', plotBoxAspectRatio)
xlabel('zscore( Log( ReactionTime ) )')
ylabel('number of samples')
dummy = [ sum(zScoreMojRT(idx,1)<-2) sum(zScoreMojRT(idx,1)>2) ];
text(-4,250, [num2str(dummy(1)) ' trials'], 'Color', 'r')
text(2.5,250, [num2str(dummy(2)) ' trials'], 'Color', 'r')
clear  idx counts dummy
        
xLimit = [min(zScoreMacRT(:,1)) , max(zScoreMacRT(:,1))] + [-1 1];
% plot MacDuff right saccade
subplot( numRows, numColumns, 2 ), hold on, box on
title('MacDuff right saccade')
idx = find(zScoreMacRT(:,2)>0);
[counts,~] = hist(zScoreMacRT(idx,1), xbins);
bar(xbins, counts, 'FaceColor', [153 255 204]/255)
plot([0 0] , yLimit, '-- k')
plot([2 2] , yLimit, '-- r', [-2 -2] , yLimit, '-- r')
set(gca, 'YLim', yLimit, 'XLim', xLimit, 'XTick', sort(unique( [floor(xLimit(1)):2:ceil(xLimit(2)) 2 -2] )), 'PlotBoxAspectRatio', plotBoxAspectRatio)
xlabel('zscore( Log( ReactionTime ) )')
ylabel('number of samples')
dummy = [ sum(zScoreMacRT(idx,1)<-2) sum(zScoreMacRT(idx,1)>2) ];
text(-4,250, [num2str(dummy(1)) ' trials'], 'Color', 'r')
text(2.5,250, [num2str(dummy(2)) ' trials'], 'Color', 'r')
clear  idx counts dummy

% plot MacDuff left saccade
subplot( numRows, numColumns, 4 ), hold on, box on
title('MacDuff left saccade')
idx = find(zScoreMacRT(:,2)<0);
[counts,~] = hist(zScoreMacRT(idx,1), xbins);
bar(xbins, counts, 'FaceColor', [153 255 204]/255)
plot([0 0] , yLimit, '-- k')
plot([2 2] , yLimit, '-- r', [-2 -2] , yLimit, '-- r')
set(gca, 'YLim', yLimit, 'XLim', xLimit, 'XTick', sort(unique( [floor(xLimit(1)):2:ceil(xLimit(2)) 2 -2] )), 'PlotBoxAspectRatio', plotBoxAspectRatio)
xlabel('zscore( Log( ReactionTime ) )')
ylabel('number of samples')
dummy = [ sum(zScoreMacRT(idx,1)<-2) sum(zScoreMacRT(idx,1)>2) ];
text(-4,250, [num2str(dummy(1)) ' trials'], 'Color', 'r')
text(2.5,250, [num2str(dummy(2)) ' trials'], 'Color', 'r')
clear  idx counts dummy

%% plot all reaction time samples for two monkeys
figure 

xLimit = [min(zScoreMojRT(:,1)) , max(zScoreMojRT(:,1))] + [-1 1];
subplot( 1, 2, 1 ), hold on, box on
title('Mojo all saccade')
[counts,~] = hist(zScoreMojRT(:,1), xbins);
bar(xbins, counts, 'FaceColor', [153 255 204]/255)
plot([0 0] , yLimit, '-- k')
plot([2 2] , yLimit, '-- r', [-2 -2] , yLimit, '-- r')
set(gca, 'YLim', yLimit, 'XLim', xLimit, 'XTick', floor(xLimit(1)):2:ceil(xLimit(2)), 'PlotBoxAspectRatio', plotBoxAspectRatio)
xlabel('zscore( Log( ReactionTime ) )')
ylabel('number of samples')
dummy = [ sum(zScoreMojRT(:,1)<-2) sum(zScoreMojRT(:,1)>2) ];
text(-4,250, [num2str(dummy(1)) ' trials'], 'Color', 'r')
text(2.5,250, [num2str(dummy(2)) ' trials'], 'Color', 'r')
clear  counts dummy

xLimit = [min(zScoreMacRT(:,1)) , max(zScoreMacRT(:,1))] + [-1 1];
subplot( 1, 2, 2 ), hold on, box on
title('MacDuff all saccade')
[counts,~] = hist(zScoreMacRT(:,1), xbins);
bar(xbins, counts, 'FaceColor', [153 255 204]/255)
plot([0 0] , yLimit, '-- k')
plot([2 2] , yLimit, '-- r', [-2 -2] , yLimit, '-- r')
set(gca, 'YLim', yLimit, 'XLim', xLimit, 'XTick', floor(xLimit(1)):2:ceil(xLimit(2)), 'PlotBoxAspectRatio', plotBoxAspectRatio)
xlabel('zscore( Log( ReactionTime ) )')
ylabel('number of samples')
dummy = [ sum(zScoreMacRT(:,1)<-2) sum(zScoreMacRT(:,1)>2) ];
text(-4,250, [num2str(dummy(1)) ' trials'], 'Color', 'r')
text(2.5,250, [num2str(dummy(2)) ' trials'], 'Color', 'r')
clear  counts dummy



