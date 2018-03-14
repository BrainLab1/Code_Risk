


clear, clc
dbstop if error 

spike_dir   = 'X:\data\Add PSTH to Aligned Spk\'; % 'Spike_PSTH_ValueTask\Spike_PSTH_CondCodeAdded_ValueTask\';

Monkey = 'Mojo';         % select among 'Mojo', 'MacDuff' 

% Collect And Store Spike Counts if it was asked
CollectAndStoreSpikeCounts = 1;

alignedEvent = '(cue)';     % select among 'cue', 'reward' and    
groupingCriteria = 'expecter_reward & RewardVariance';  % see all possible conditions under GroupTrials()
timeInterval = [0,0.2];
minRepetitionPerCondition = 5;
ZScoreEachTuningCurve = 1;

color_code = hsv(9);   % ***** use specific color for each task condition; 36 colors

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

if CollectAndStoreSpikeCounts

idx = 1;
% read out the session files one by one
while idx <= length(allSessions)
    
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
    
    save(['X:\Report 11\Spike Counts\' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']_' ...
           groupingCriteria '_minRepetition' num2str(minRepetitionPerCondition) '\' allSessions(idx).name], 'spkDataArr1', 'spkDataArr2')
    clear spkDataArr* decision unitNames*
    idx = idx+1;
end

end % end of CollectAndStoreSpikeCounts

%%
SpkCountsPath = ['X:\Report 11\Spike Counts\' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']_' groupingCriteria '_minRepetition' num2str(minRepetitionPerCondition) '\'];

% read out list of spike count files
sesListSpkCnt = dir(SpkCountsPath);
% remove '.' and '..'
sesListSpkCnt = sesListSpkCnt(3:end);

% select the right data files based on above options
idx = 1;
ind = [];
while idx <= length(sesListSpkCnt)
    % find the file names which do not match the settings above
    if ( isempty(strfind(sesListSpkCnt(idx).name,Monkey)) || ...
         isempty(strfind(sesListSpkCnt(idx).name,alignedEvent)))
       ind = [ind,idx];
    end
    idx = idx+1;
end
% remove the unmatched file names/data files/sessions
sesListSpkCnt(ind) = [];
clear idx ind

ses = 1;
while ses <= length(sesListSpkCnt)
    ses
    % load this sesion file; this will load spkDataArr1 and spkDataArr2
    load([SpkCountsPath sesListSpkCnt(ses).name])
    
    % collect mean spike count of all the units, separately for each group 
    if ~exist('allUnitsSpkCountArr1','var')
        allUnitsSpkCountArr1 = []; % each row will be the spike count for one unit
    end
    if ~exist('allUnitsSpkCountArr2','var')
        allUnitsSpkCountArr2 = []; % each row will be the spike count for one unit
    end
    
    % collect spike counts from all the units on array 1
    unNames = fieldnames(spkDataArr1);
    for un=1:length(unNames)
        if ~isempty (spkDataArr1.(unNames{un}))
            switch ZScoreEachTuningCurve
                case 1
                    allUnitsSpkCountArr1 = [allUnitsSpkCountArr1; zscore( cellfun(@(x) nanmean(x), {spkDataArr1.(unNames{un}).SpikeCount}) ) ]; % each row will be the spike count for one unit
                case 0
                    allUnitsSpkCountArr1 = [allUnitsSpkCountArr1; ( cellfun(@(x) nanmean(x), {spkDataArr1.(unNames{un}).SpikeCount}) ) ]; % each row will be the spike count for one unit
            end
        end
    end
    clear unNames un
    
    % collect spike counts from all the units on array 2
    unNames = fieldnames(spkDataArr2);
    for un=1:length(unNames)
        if ~isempty (spkDataArr2.(unNames{un}))
            switch ZScoreEachTuningCurve
                case 1
                    allUnitsSpkCountArr2 = [allUnitsSpkCountArr2; zscore( cellfun(@(x) nanmean(x), {spkDataArr2.(unNames{un}).SpikeCount}) ) ]; % each row will be the spike count for one unit
                case 0
                    allUnitsSpkCountArr2 = [allUnitsSpkCountArr2; ( cellfun(@(x) nanmean(x), {spkDataArr2.(unNames{un}).SpikeCount}) ) ]; % each row will be the spike count for one unit
            end
        end
    end
    clear unNames un
    
    clear spkDataArr*
    ses = ses + 1;
end
clear ses

%% plot the results

switch groupingCriteria
    case 'Variance'
        figure('Name', Monkey)
        subplot(2,1,1), box on, title('dlPFC')
        plot3(allUnitsSpkCountArr1(:,1) , allUnitsSpkCountArr1(:,2) , allUnitsSpkCountArr1(:,3),'.')
        xlabel('Var = 3') 
        ylabel('Var = 6') 
        zlabel('Var = 9') 

        subplot(2,1,2), box on, title('7a')
        plot3(allUnitsSpkCountArr2(:,1) , allUnitsSpkCountArr2(:,2) , allUnitsSpkCountArr2(:,3),'.')
        xlabel('Var = 3') 
        ylabel('Var = 6') 
        zlabel('Var = 9') 

        
    case 'expecter_reward & RewardVariance'
        fg = figure('Name', Monkey);
        paraVal = [3 0; 3 1; 3 4; 6 0; 6 1; 6 4; 9 0; 9 1; 9 4];
        colorCode = [0 1 1 
                     0 0 1
                     0 0 0];
       
       % ----- plot array 1 
        subplot(2,1,1), title(['dlPFC ' alignedEvent ' [' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']']), hold on, box on
        xlabel('EV'); 
        switch ZScoreEachTuningCurve
            case 1
                ylabel('grand average Z-Score(spike count) (+/- sem)')
            case 0
                ylabel('grand average spike count (+/- sem)')
        end
        txt = ['N = ' num2str(size(allUnitsSpkCountArr1,1))];
        for vr = [0 1 4]
            % initialize cuvVarX
            eval(['crvVar' num2str(vr) ' = [];']);
            eval(['semVar' num2str(vr) ' = [];']);
            % find the group parameters with this variance value
            idxVar = find(paraVal(:,2) == vr);
            for ev = 1:length(idxVar) % these are Ev 3, 6, 9
                    eval([ 'crvVar' num2str(vr) ' = [crvVar' num2str(vr) ', '  num2str(mean(allUnitsSpkCountArr1(:,idxVar(ev)))) ']' ]);
                    eval([ 'semVar' num2str(vr) ' = [semVar' num2str(vr) ', '  num2str(std(allUnitsSpkCountArr1(:,idxVar(ev)))/size(allUnitsSpkCountArr1,1)) ']' ]);    
            end
            clear idxVar ev
            eval([ 'errorbar( 1:3, crvVar' num2str(vr) ', semVar' num2str(vr) ', ''Color'', colorCode(sqrt(vr)+1,:) )'])
            txt = [txt; {['\color[rgb]{' num2str(colorCode(sqrt(vr)+1, :)) '} Var = ' num2str(vr)]}];
        end
        set(gca, 'XTick', 1:3, 'XTickLabel', [{'3'},{'6'},{'9'}], 'XLim', [0 4])
        text(2,mean(ylim),txt)
        clear sem* crv*
        aa = allUnitsSpkCountArr1';
        [p,tbl,stats] = anovan(aa(:), repmat(paraVal, size(allUnitsSpkCountArr1,1), 1), 'model', 'interaction', 'varnames', {'EV','Var'}, 'display', 'off');
        clear aa
        txtANOVA = [{['pEV = ' num2str(p(1))]} , {['pVar = ' num2str(p(2))]} , {['pEV*Var = ' num2str(p(3))]}];
        text(3,0,txtANOVA)
        clear txtANOVA p tbl stats
        
       % ----- plot array 2 
        subplot(2,1,2), title(['7a ' alignedEvent ' [' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']']), hold on, box on
        xlabel('EV'); 
        switch ZScoreEachTuningCurve
            case 1
                ylabel('grand average Z-Score(spike count) (+/- sem)')
            case 0
                ylabel('grand average spike count (+/- sem)')
        end
        txt = ['N = ' num2str(size(allUnitsSpkCountArr2,1))];
        for vr = [0 1 4]
            % initialize cuvVarX
            eval(['crvVar' num2str(vr) ' = [];']);
            eval(['semVar' num2str(vr) ' = [];']);
            % find the group parameters with this variance value
            idxVar = find(paraVal(:,2) == vr);
            for ev = 1:length(idxVar) % these are Ev 3, 6, 9
                    eval([ 'crvVar' num2str(vr) ' = [crvVar' num2str(vr) ', '  num2str(mean(allUnitsSpkCountArr2(:,idxVar(ev)))) ']' ]);
                    eval([ 'semVar' num2str(vr) ' = [semVar' num2str(vr) ', '  num2str(std(allUnitsSpkCountArr2(:,idxVar(ev)))/size(allUnitsSpkCountArr2,1)) ']' ]);    
            end
            clear idxVar ev
            eval([ 'errorbar( 1:3, crvVar' num2str(vr) ', semVar' num2str(vr) ', ''Color'', colorCode(sqrt(vr)+1,:) )'])
            txt = [txt; {['\color[rgb]{' num2str(colorCode(sqrt(vr)+1, :)) '} Var = ' num2str(vr)]}];
        end
        set(gca, 'XTick', 1:3, 'XTickLabel', [{'3'},{'6'},{'9'}], 'XLim', [0 4])
        text(2,mean(ylim),txt)
        clear sem* crv*        
        aa = allUnitsSpkCountArr2';
        [p,tbl,stats] = anovan(aa(:), repmat(paraVal, size(allUnitsSpkCountArr2,1), 1), 'model', 'interaction', 'varnames', {'EV','Var'}, 'display', 'off');
        clear aa
        txtANOVA = [{['pEV = ' num2str(p(1))]} , {['pVar = ' num2str(p(2))]} , {['pEV*Var = ' num2str(p(3))]}];
        text(3,0,txtANOVA)
        clear txtANOVA p tbl stats
        
        flName = [Monkey '_MeanPopFiring_' regexprep(groupingCriteria,' & ','-')];
        saveas(fg, [SpkCountsPath '\Results\' flName], 'fig')
        saveas(fg, [SpkCountsPath '\Results\' flName], 'bmp')
        close(fg); clear fg flName
end            





