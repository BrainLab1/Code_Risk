% This scripts collects spike counts for units and stores them on the
% hard drive.

% Notice that to run this script, you need the spike count data to be
% already stored on the hard drive. You should use PlotTuningCurve.m to get
% the spike count files.

clear, clc
dbstop if error 

spike_dir   = 'X:\data\Add PSTH to Aligned Spk\'; % 'Spike_PSTH_ValueTask\Spike_PSTH_CondCodeAdded_ValueTask\';

Monkey = 'Mojo';         % select among 'Mojo', 'MacDuff' 

alignedEvent = '(cue)';     % select among 'cue', 'reward' and    
groupingCriteria = 'expecter_reward & RewardVariance';  % see all possible conditions under GroupTrials()
timeInterval = [0,0.8];
minRepetitionPerCondition = -5;

%% read out list of all the available sesions
allSessions = dir(spike_dir);
% remove the '.' and '..' from the list of sessions
allSessions = allSessions(3:end);  

%% select the right data files based on above options
idx = 1;
ind = [];
while idx <= length(allSessions)
    % find the file names which do not match the settings above
    if ( isempty(strfind(allSessions(idx).name,Monkey)) || ...
         isempty(strfind(allSessions(idx).name,alignedEvent)))
       ind = [ind,idx];
    end
    idx = idx+1;
end
% remove the unmatched file names/data files/sessions
allSessions(ind) = [];
clear idx ind

%% Collect And Store Spike Counts if it was asked
   
idx = 1;
% read out the session files one by one
while idx <= length(allSessions)
    
    display(['ses  ' num2str(idx)])
    % load the spikeTrials for current sessions
    load ([spike_dir allSessions(idx).name]);
    
    % get the list of the unit names for the current session
    unitNames = fieldnames(spikeTrials.PSTH);
    validFields = cellfun( @(x) strcmp(x(1:4),'elec') , unitNames);
    unitNames = unitNames(validFields);
    clear validFields
    
    % split units by array number 
    array1Idx = find( cellfun( @(x) strcmp(x(1:5),'elec1') , unitNames) );
    array2Idx = find( cellfun( @(x) strcmp(x(1:5),'elec2') , unitNames) );
    unitNamesArr1 = unitNames(array1Idx); 
    unitNamesArr2 = unitNames(array2Idx); 
    clear array1Idx array2Idx
    
    % get the spike count for units on array1
    for un = 1:length(unitNamesArr1)
        display(['[ses unArray1] = [' num2str(idx) ' , ' num2str(un) ']'])
        % group trials based on the groupingCriteria
        [output] = GroupTrials(spikeTrials.cfg.event, groupingCriteria);
        
        % find for which trials in each group, the unit was active 
        output = FindRepetitionPerConditionForUnit (groupingCriteria, spikeTrials.cfg.event, spikeTrials.PSTH.(unitNamesArr1{un}));

        % check if the unit has enough repetition per condition/group
        [output, decision] = SelectOrRejectUnit( output, 'RepetitionPerCondition', struct('MinRepPerCnd',minRepetitionPerCondition) );
        
        switch decision
            case 1
                % add spike counts for trials of each group
                [output] = GetSpikeCount (output, spikeTrials.PSTH.(unitNamesArr1{un}), 'cue', timeInterval);
                spkDataArr1.(unitNamesArr1{un}) = output;
            case 0
                 spkDataArr1.(unitNamesArr1{un}) = [];
        end
                
        clear output 
    end
    clear un
    
    
    % get the spike count for units on array2
    for un = 1:length(unitNamesArr2)
        display(['[ses unArray2] = [' num2str(idx) ' , ' num2str(un) ']'])
        % group trials based on the groupingCriteria
        [output] = GroupTrials(spikeTrials.cfg.event, groupingCriteria);
        
        % find for which trials in each group, the unit was active 
        output = FindRepetitionPerConditionForUnit (groupingCriteria, spikeTrials.cfg.event, spikeTrials.PSTH.(unitNamesArr2{un}));

        % check if the unit has enough repetition per condition/group
        [output, decision] = SelectOrRejectUnit( output, 'RepetitionPerCondition', struct('MinRepPerCnd',minRepetitionPerCondition) );
        
        switch decision
            case 1
                % add spike counts for trials of each group
                [output] = GetSpikeCount (output, spikeTrials.PSTH.(unitNamesArr2{un}), 'cue', timeInterval);
                spkDataArr2.(unitNamesArr2{un}) = output;
            case 0
                 spkDataArr2.(unitNamesArr2{un}) = [];
        end
                
        clear output 
    end
    clear un
    
    save(['X:\Test\Spike Counts\' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']_' ...
           groupingCriteria '_minRepetition' num2str(minRepetitionPerCondition) '\' allSessions(idx).name], 'spkDataArr1', 'spkDataArr2')
    clear spkDataArr* decision unitNames*
    idx = idx+1;
end


