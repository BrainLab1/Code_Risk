
% this script plots and stores tuning curve of GLM-significant active units for both areas and all GLM factors to 9 [EV Var] 
% task conditions.
% To run this script we need to know which of the active units had significant GLM cofficients. To get this you need to 
% run: 1- CollectAndStoreSpikeCounts and 2- R13_SaveSigUnits_ScatterPlotGLMCoeff.m


clear, clc
dbstop if error 

Monkey = 'Mojo';         % select among 'Mojo', 'MacDuff' 

alignedEvent = '(cue)';     % select among 'cue', 'reward' and    
groupingCriteria = 'expecter_reward & RewardVariance'; %'CuePos & CurrEV & CurrVar';  % see all possible conditions under GroupTrials()
timeInterval = [0,0.4];
minRepetitionPerCondition = -5;
ZScoreEachTuningCurve = 0;

colorCode = [[0 1 1]; ...
             [0 0 1]; ...
             [0 0 0]];

SpkCountsPath = ['X:\Test\Spike Counts\' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']_' groupingCriteria '_minRepetition' num2str(minRepetitionPerCondition) '\'];
% % % % ParaNames = {'ev', 'var', 'int'};
Area = {'dlPFC','7a'};

numRows = 2;
numColumns = 3;

%% select the right data files based on above options
spike_dir   = 'X:\data\Add PSTH to Aligned Spk\'; %'D:\shared\data\Add PSTH to Aligned Spk\'; 

% read out list of all the available sesions
allSessions = dir(spike_dir);
% remove the '.' and '..' from the list of sessions
allSessions = allSessions(3:end);  

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

%%
idx = 11;

% read out the session files one by one  
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

% % % %         % read out the event structure for this unit active trials
% % % %         thisUnitEven = spikeTrials.cfg.event(spikeTrials.PSTH.(unitNamesArr1{un}).active_eventStrctTrls);
% % % %         % group trials based on the groupingCriteria
% % % %         [output] = GroupTrials(thisUnitEven, groupingCriteria);

        % group trials based on the groupingCriteria
        [output] = GroupTrials(spikeTrials.cfg.event, groupingCriteria);
        
        % find for which trials in each group, the unit was active 
        output = FindRepetitionPerConditionForUnit (groupingCriteria, spikeTrials.cfg.event, spikeTrials.PSTH.(unitNamesArr1{un}));

        % add repetition per condition to the output and check if the unit has enough repetition per condition/group
        [output, ~] = SelectOrRejectUnit( output, 'RepetitionPerCondition', struct('MinRepPerCnd',minRepetitionPerCondition) );
         
        % add spike counts for trials of each group
        [output] = GetSpikeCount (output, spikeTrials.PSTH.(unitNamesArr1{un}), 'cue', timeInterval);
        spkDataArr1.(unitNamesArr1{un}) = output;
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

        % add repetition per condition to the output and check if the unit has enough repetition per condition/group
        [output, ~] = SelectOrRejectUnit( output, 'RepetitionPerCondition', struct('MinRepPerCnd',minRepetitionPerCondition) );

        % add spike counts for trials of each group
        [output] = GetSpikeCount (output, spikeTrials.PSTH.(unitNamesArr2{un}), 'cue', timeInterval);
        spkDataArr2.(unitNamesArr2{un}) = output;
                
        clear output 
    end
    clear un
    
%%
% initialize figure
fgCount = 1;
fg = figure('Name',[Monkey '_' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']'], 'Position', [40 70 0 0] + 0.8*get(0,'ScreenSize'));
subPltCount = 0;

% plot array1 tuning curves -----------------------------
for u = 1:length(unitNamesArr1)
    subPltCount = subPltCount+1;
    % read out unit's spike counts, already grouped
    output = spkDataArr1.(unitNamesArr1{u});
    
    % get the parameter values for each group
    paraVal = cell2mat({output.Value}');

    % find the right place to plot the tuning curve!
    if subPltCount > (numRows*numColumns)
        saveas(fg,[SpkCountsPath 'Results\' Monkey '_dlPFC_' num2str(fgCount)], 'bmp')
        close(fg)
        clear fg
        fg = figure('Name',[Monkey '_' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']'], 'Position', [40 70 0 0] + 0.8*get(0,'ScreenSize'));
        fgCount = fgCount+1;
        subPltCount = 1;
    end            
    figure(fg), 
    subplot(numRows, numColumns, subPltCount), hold on, box on
    title(regexprep(unitNamesArr1{u},'_','-'))
    xlabel('EV', 'FontWeight','Bold'); 
    ylabel('Mean Firing Rate', 'FontWeight','Bold')
    
    txt = [];
    for vr = [0 1 4]
        % initialize cuvVarX
        eval(['crvVar' num2str(vr) ' = [];']);
        eval(['semVar' num2str(vr) ' = [];']);
        % find the group parameters with this variance value
        idxVar = find(paraVal(:,2) == vr);
        for eevv = 1:length(idxVar) % these are Ev 3, 6, 9
            eval([ 'crvVar' num2str(vr) ' = [crvVar' num2str(vr) ', '  num2str(mean(output(idxVar(eevv)).SpikeCount)) '];' ]);
            eval([ 'semVar' num2str(vr) ' = [semVar' num2str(vr) ', '  num2str(std(output(idxVar(eevv)).SpikeCount)/length(output(idxVar(eevv)).SpikeCount)) '];' ]);    
        end
        clear idxVar eevv
        eval([ 'errorbar( 1:3, crvVar' num2str(vr) ', semVar' num2str(vr) ', ''Color'', colorCode(sqrt(vr)+1,:) )'])
        txt = [txt; {['\color[rgb]{' num2str(colorCode(sqrt(vr)+1, :)) '} Var = ' num2str(vr)]}];
    end
    set(gca, 'XTick', 1:3, 'XTickLabel', [{'3'},{'6'},{'9'}], 'XLim', [0 4])
    text(2,mean(ylim),txt)
    clear sem* crv* paraVal txt vr
end
clear u
saveas(fg,[SpkCountsPath 'Results\' Monkey '_dlPFC_' num2str(fgCount)], 'bmp')
close(fg)
clear fg fgCount
    
% plot array2 tuning curves -----------------------
% initialize figure
fgCount = 1;
fg = figure('Name',[Monkey '_' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']'], 'Position', [40 70 0 0] + 0.8*get(0,'ScreenSize'));
subPltCount = 0;

for u = 1:length(unitNamesArr2)
    subPltCount = subPltCount+1;
    % read out unit's spike counts, already grouped
    output = spkDataArr2.(unitNamesArr2{u});
    
    % get the parameter values for each group
    paraVal = cell2mat({output.Value}');

    % find the right place to plot the tuning curve!
    if subPltCount > (numRows*numColumns)
        saveas(fg,[SpkCountsPath 'Results\' Monkey '_7a_' num2str(fgCount)], 'bmp')
        close(fg)
        clear fg
        fg = figure('Name',[Monkey '_' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']'], 'Position', [40 70 0 0] + 0.8*get(0,'ScreenSize'));
        fgCount = fgCount+1;
        subPltCount = 1;
    end            
    figure(fg), 
    subplot(numRows, numColumns, subPltCount), hold on, box on
    title(regexprep(unitNamesArr2{u},'_','-'))
    xlabel('EV', 'FontWeight','Bold'); 
    ylabel('Mean Firing Rate', 'FontWeight','Bold')
    
    txt = [];
    for vr = [0 1 4]
        % initialize cuvVarX
        eval(['crvVar' num2str(vr) ' = [];']);
        eval(['semVar' num2str(vr) ' = [];']);
        % find the group parameters with this variance value
        idxVar = find(paraVal(:,2) == vr);
        for eevv = 1:length(idxVar) % these are Ev 3, 6, 9
            eval([ 'crvVar' num2str(vr) ' = [crvVar' num2str(vr) ', '  num2str(mean(output(idxVar(eevv)).SpikeCount)) '];' ]);
            eval([ 'semVar' num2str(vr) ' = [semVar' num2str(vr) ', '  num2str(std(output(idxVar(eevv)).SpikeCount)/length(output(idxVar(eevv)).SpikeCount)) '];' ]);    
        end
        clear idxVar eevv
        eval([ 'errorbar( 1:3, crvVar' num2str(vr) ', semVar' num2str(vr) ', ''Color'', colorCode(sqrt(vr)+1,:) )'])
        txt = [txt; {['\color[rgb]{' num2str(colorCode(sqrt(vr)+1, :)) '} Var = ' num2str(vr)]}];
    end
    set(gca, 'XTick', 1:3, 'XTickLabel', [{'3'},{'6'},{'9'}], 'XLim', [0 4])
    text(2,mean(ylim),txt)
    clear sem* crv* paraVal txt vr
end
clear u
saveas(fg,[SpkCountsPath 'Results\' Monkey '_7a_' num2str(fgCount)], 'bmp')
close(fg)
clear fg
    
        