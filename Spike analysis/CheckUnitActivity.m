

% this script checks compares units have ave rage firing rate within the StabilityCheckTimeInterval to
% AveSpkCountThreshold, and split them into ActiveSesUnitsArrX and SilentSesUnitsArrX, stores list of 
% active and silent units and session names on the hard drive

clear, clc
dbstop if error 

spike_dir   = 'X:\data\Add PSTH to Aligned Spk\'; % 'Spike_PSTH_ValueTask\Spike_PSTH_CondCodeAdded_ValueTask\';

Monkey = 'MacDuff';         % select among 'Mojo', 'MacDuff' 

alignedEvent = '(cue)';     % select among 'cue', 'reward' and    
groupingCriteria = 'CuePos & CurrEV & CurrVar';  % see all possible conditions under GroupTrials()
timeInterval = [0,0.2];
minRepetitionPerCondition = 5;
ZScoreEachTuningCurve = 1;

color_code = hsv(9);   % ***** use specific color for each task condition; 36 colors

ApplyGLMonPopulation = 0;
ApplyGLMonIndividualUnits = 1 ;

StabilityCheckTimeInterval = [-0.5 1];
AveSpkCountThreshold = 5;

ActiveSesUnitsArr1 = [];
SilentSesUnitsArr1 = [];
ActiveSesUnitsArr2 = [];
SilentSesUnitsArr2 = [];
%% read out list of all the available sesions
allSessions = dir(spike_dir);
% remove the '.' and '..' from the list of sessions
allSessions = allSessions(3:end);  
% select the right data files based on above options
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
        aa = spikeTrials.PSTH.(unitNamesArr1{un}).SpikeTimeByTrial;
        % in each trial, find spikes which are happening within the StabilityCheckTimeInterval
        bb = cellfun(@(x)  x( find((x>=StabilityCheckTimeInterval(1)) .* (x<=StabilityCheckTimeInterval(2))) ) , aa , 'UniformOutput' , 0);
        % in each trial, get the number of spikes within the StabilityCheckTimeInterval
        cc = cellfun(@(x) length(x), bb);
        
        switch mean(cc)>=AveSpkCountThreshold
            case 1
                ActiveSesUnitsArr1  = [ActiveSesUnitsArr1; [{allSessions(idx).name} , {unitNamesArr1(un)} , {mean(cc)}]];
            case 0
                SilentSesUnitsArr1  = [SilentSesUnitsArr1; [{allSessions(idx).name} , {unitNamesArr1(un)} , {mean(cc)}]];
        end
        clear aa bb cc    
    end
    clear un
    
    
    % get the spike count for units on array2
    for un = 1:length(unitNamesArr2)
        display(['[ses unArray2] = [' num2str(idx) ' , ' num2str(un) ']'])
        aa = spikeTrials.PSTH.(unitNamesArr2{un}).SpikeTimeByTrial;
        % in each trial, find spikes which are happening within the StabilityCheckTimeInterval
        bb = cellfun(@(x)  x( find((x>=StabilityCheckTimeInterval(1)) .* (x<=StabilityCheckTimeInterval(2))) ) , aa , 'UniformOutput' , 0);
        % in each trial, get the number of spikes within the StabilityCheckTimeInterval
        cc = cellfun(@(x) length(x), bb);
        
        switch mean(cc)>=AveSpkCountThreshold
            case 1
                ActiveSesUnitsArr2  = [ActiveSesUnitsArr2; [{allSessions(idx).name} , {unitNamesArr2(un)} , {mean(cc)}]];
            case 0
                SilentSesUnitsArr2  = [SilentSesUnitsArr2; [{allSessions(idx).name} , {unitNamesArr2(un)} , {mean(cc)}]];
        end
        clear aa bb cc    
    end
    clear un
    
    clear spikeTrials
    idx = idx+1;
end


save(['X:\data\' Monkey '_UnitStability_' alignedEvent '[' regexprep(num2str(StabilityCheckTimeInterval(1)),'0.','0,') ' ' regexprep(num2str(StabilityCheckTimeInterval(2)),'0.','0,') ']_AveSpkCntThreshold' num2str(AveSpkCountThreshold)], ...
     'ActiveSesUnitsArr1','ActiveSesUnitsArr2','SilentSesUnitsArr1','SilentSesUnitsArr2','AveSpkCountThreshold','StabilityCheckTimeInterval','Monkey')


