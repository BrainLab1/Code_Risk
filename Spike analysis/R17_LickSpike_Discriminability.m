
% initiated on 12.05.2018

clear
clc
dbstop if error

% General settings
Monkey = 'Mojo';  % 'Moj' or 'Mac'

% for behavioral data files
bhv_data_folder             = 'Z:\data\Risk\Behavior\';
lickGroupingCriteria       = 'expecter_reward & RewardVariance' ; 
numConditions = 9; 
lickFs = 1000;   % assume 1K sampling rate
switch Monkey
    case 'MacDuff'
        eventTimeInterval = {[0.75 0.95]};  % for lick signal
        alignEvent = {'cue'};  
    case 'Mojo'
        eventTimeInterval = {[0.4 0.8]};
        alignEvent = {'cue'};  
end

% for spike data files
spike_dir   = 'X:\data\Add PSTH to Aligned Spk\'; % 'Spike_PSTH_ValueTask\Spike_PSTH_CondCodeAdded_ValueTask\';
spikeGroupingCriteria = 'CuePos & CurrEV & CurrVar';  % see all possible conditions under GroupTrials()
spikeAlignedEvent = '(cue)';
spikeTimeInterval = [0.2,0.4];
spikeMinRepetitionPerCondition = 5;
ZScoreEachTuningCurve = 1;
StabilityCheckTimeInterval = [-0.5 1];
AveSpkCountThreshold = 5;
         
save_dir = 'X:\Report 17\Spike vs Lick discriminability\';
DiscriminabilityType = 'Cohen''s d';
%% Read out list of all the files related to this session
allFiles = dir(bhv_data_folder);    

% remove the '.' and '..' from the folders list and filter sessions for monkey names
idx = [];
for ses = 1:numel(allFiles) % for each session
    if (numel(allFiles(ses).name) < 3) || ~(strcmp(allFiles(ses).name(5:7), Monkey(1:3)))  
        idx = [idx; ses];
    end
end
allFiles(idx) = [];
clear ses idx
%% Collect lick signal
allSesGroupedData = {};
for ses = 1:numel(allFiles) % for each session
    display(['Session ' num2str(ses)])

    % this will load new_cfg variable
    load ([bhv_data_folder allFiles(ses).name])
    
    % convert the event tructure to table
    eventTable = struct2table(new_cfg.event);
    
    % ........remove outliers....
    % find outlier trials for duration of cue presentation & reward time and reaction time
    outTrials = find(eventTable.OutCueOnOffZ2 + eventTable.OutAqToRwdZ2 + eventTable.OutRTZ2);
    clear eventTable
    % remove oulier trials from eventTable and new_cfg.event 
    new_cfg.event(outTrials) = [];
    clear outTrials
 
    % group trials for 27 task conditions
    [output] = GroupTrials(new_cfg.event, lickGroupingCriteria);   
    
    % if doesn't exist, make a new variable with the same size as the number of groups
    if ~exist('allSesSuccessLick','var')
        allSesSuccessLick = cell(length(output), length(alignEvent));
    end
    
    % convert the outlier-removed event tructure to table
    eventTable = struct2table(new_cfg.event);
    

    for tim = 1:length(alignEvent) % for each time interval
        for gr = 1:length(output)
            % get the trial index for successful trials per condition
            output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
            % get the reward variance for output(gr).SuccessTrlIdx trials
            output(gr).CurrVar = eventTable.RewardVariance(output(gr).SuccessTrlIdx);
            % get the expected reward for output(gr).SuccessTrlIdx trials
            output(gr).CurrEV = eventTable.expected_reward(output(gr).SuccessTrlIdx);

            % for successful trials of this group, collect the lick signal in the asked time interval 
            temp = [];
            cutInterval = [];
            for tr = 1:length(output(gr).SuccessTrlIdx)
                trIdx = output(gr).SuccessTrlIdx(tr);
                switch alignEvent{tim}
                    case 'cue'
                        % cut the licking signall based on DiodeEventTimes
                        cutInterval(tr,:) = eventTable.DiodeEventTime{trIdx}.CueOnset + (eventTimeInterval{tim}*lickFs);
                        tmpSig = eventTable.LickSignal{trIdx} > 2.5;
                        temp = [temp; tmpSig(cutInterval(tr,1) : cutInterval(tr,2))'];
                        clear tmpSig
                    case 'reward'
                        cutInterval(tr,:) = eventTable.RewardOnTime{trIdx}(1) + (eventTimeInterval{tim}*lickFs);
                        tmpSig = eventTable.LickSignal{trIdx} > 2.5;
                        temp = [temp; tmpSig(cutInterval(tr,1) : cutInterval(tr,2))'];
                        clear tmpSig
                end
                clear trIdx
            end
            % get the mean lick signal for each trial
            eval(['output(gr).SuccessLickInterval' num2str(tim) ' = mean(temp,2);'])
            eval(['output(gr).OrigCutSampleInterval' num2str(tim) ' = cutInterval;'])
            clear temp tr cutInterval
        end
        clear gr
        
        % ---- zscore mean lick signals in this session
        % get number of trials in each group
        eval(['sz = cellfun(@(x) length(x), {output.SuccessLickInterval' num2str(tim) '}'');'])
        % zscore mean lick across trials
        bb = struct2table(output); 
        eval(['bb.SuccessLickInterval' num2str(tim) ' = mat2cell(zscore(cell2mat({output.SuccessLickInterval' num2str(tim) '}'')), sz, 1);'])
        clear output
        output = table2struct(bb);
        clear bb sz
        
        for gr = 1:length(output)
            % collect lick signals of this group in this session's variable
            eval(['allSesSuccessLick(gr,tim) = {[allSesSuccessLick{gr,tim}; output(gr).SuccessLickInterval' num2str(tim) ']};'])
        end
        clear gr
    end
    allSesGroupedData{ses} = output;
    clear gr output eventTable new_cfg tmp_cfg temp tim
end
clear ses

%%  Get the discriminability index across conditions for the lick signal
lickDiscrMatrix = cell(size(allSesSuccessLick,2)); % length(lickDiscrMatrix) = number of time intervals
for tim = 1:length(lickDiscrMatrix)
   % initialize the discriminablity matrix
    lickDiscrMatrix{tim} = nan(size(allSesSuccessLick,1)); 
% % % %     lickDiscrMatrix{tim} = cell2mat( cellfun(@(x,y)   cell2mat(cellfun(@(x1,y)  (mean(x1)-mean(y))/sqrt((var(x1)/length(x1))+(var(y)/length(y)))  , repmat({x},size(allSesSuccessLick,1),1),  allSesSuccessLick(:,tim), 'UniformOutput', 0))'  , allSesSuccessLick(:,tim), 'UniformOutput', 0) );
    switch DiscriminabilityType 
        case 'Cohen''s d'
            lickDiscrMatrix{tim} = cell2mat( cellfun(@(x,y)   cell2mat(cellfun(@(x1,y)  (mean(x1)-mean(y))/sqrt( 0.5*((nanstd(x1)^2)+(nanstd(y)^2)) )  , repmat({x},size(allSesSuccessLick,1),1),  allSesSuccessLick(:,tim), 'UniformOutput', 0))'  , allSesSuccessLick(:,tim), 'UniformOutput', 0) );
    end
end
clear tim

%% Spike data
Array_pos = {'Anterior','Posterior'};% Array1 is anterior and Array2 is posterior
eval(['[newArray1, newArray2] = new_electrodepinout' Monkey '(''chan'');']) 

% here we use spike count data which is already stored on the hard drive
SpkCountsPath = ['X:\Report 13\Spike Counts\' spikeAlignedEvent(2:end-1) '[' num2str(spikeTimeInterval(1)) ',' num2str(spikeTimeInterval(2)) ']_' spikeGroupingCriteria '_minRepetition' num2str(spikeMinRepetitionPerCondition) '\'];

% load list of active units for this monkey
load(['X:\data\' Monkey '_UnitStability_' spikeAlignedEvent '[' regexprep(num2str(StabilityCheckTimeInterval(1)),'0.','0,') ' ' regexprep(num2str(StabilityCheckTimeInterval(2)),'0.','0,') ']_AveSpkCntThreshold' num2str(AveSpkCountThreshold)],...
      'ActiveSesUnitsArr1', 'ActiveSesUnitsArr2');

sesListSpkCnt = unique([ActiveSesUnitsArr1(:,1);ActiveSesUnitsArr2(:,1)]);

ses = 1;
popUnitIdxArr1 = 0;
popUnitIdxArr2 = 0;
while ses <= length(sesListSpkCnt)
    ses
    % get the list of active units in this session on array 1 and 2
    thisSesActiveUnits = [ ActiveSesUnitsArr1( find( strcmp(ActiveSesUnitsArr1(:,1),sesListSpkCnt{ses}) ) , 2); ...
                           ActiveSesUnitsArr2( find( strcmp(ActiveSesUnitsArr2(:,1),sesListSpkCnt{ses}) ) , 2)];
    
    % load this sesion file; this will load spkDataArr1 and spkDataArr2
    load([SpkCountsPath sesListSpkCnt{ses}])
    % merge spkDataArr1 and spkDataArr2 which were made using old electrode-channel assignment  
    spkData = [];
    fnames = fieldnames(spkDataArr1);
    for i=1:length(fnames)
        spkData.(fnames{i}) = spkDataArr1.(fnames{i});
    end
    clear i fames spkDataArr1
    fnames = fieldnames(spkDataArr2);
    for i=1:length(fnames)
        spkData.(fnames{i}) = spkDataArr2.(fnames{i});
    end
    clear i fames spkDataArr2

    % remove empty structures (it happens for some units in Mojo's data!!)
    unitNames = fieldnames(spkData);
    reovemUnits = find(structfun(@isempty, spkData));
    if ~isempty(reovemUnits)
        spkData = rmfield(spkData,unitNames(reovemUnits));
        for iii=1:length(reovemUnits)
             thisSesActiveUnits( find(cellfun(@(x) strcmp(x,unitNames{reovemUnits(iii)}), thisSesActiveUnits)) ) = [];
        end
        clear reovemUnits unitNames iii
    end
    
    % collect mean spike count of all the units, separately for each group 
    if ~exist('allUnitTrlsSpkCountArr1','var')
        allUnitTrlsSpkCountArr1 = []; % each row will be the spike count for one unit
        allUnitTrlsParaValArr1 = [];
        allUnitTrlsUnitIdxArr1 = [];
        allUnitGroupsArr1 = [];
    end
    if ~exist('allUnitTrlsSpkCountArr2','var')
        allUnitTrlsSpkCountArr2 = []; % each row will be the spike count for one unit
        allUnitTrlsParaValArr2 = [];
        allUnitTrlsUnitIdxArr2 = [];
        allUnitGroupsArr2 = [];
    end
    
    % collect spike counts from all the units 
    for un=1:length(thisSesActiveUnits)
        % find new array assignment for this channel
        charIdx = cell2mat(strfind(thisSesActiveUnits{un},'ch_'));
        chanNum = str2num(thisSesActiveUnits{un}{1,1}(charIdx+3:end));
        if ismember(chanNum, newArray1)
            thisActiveChanArrayIdx = 1;
        elseif ismember(chanNum, newArray2)
            thisActiveChanArrayIdx = 2;
        else
            display(['Erroe: Channel ' num2str(chanNum) ' was not assigned to either arrays!!!!!'])    
        end
        clear charIdx chanNum 
        
        % zscore spike counts
        tmp = cell2mat({spkData.(thisSesActiveUnits{un,1}{1,1}).SpikeCount}');
        aa = mat2cell ( (tmp - nanmean(tmp))/nanstd(tmp) , cellfun(@(x) length(x), {spkData.(thisSesActiveUnits{un,1}{1,1}).SpikeCount}'), 1);
        clear tmp
        
        % merge trial groups for lwft and right cue location
        bb =  cell2mat({spkData.(thisSesActiveUnits{un,1}{1,1}).Value}');
        uniqEVandVar = unique(bb(:,2:3),'rows');
        tempOutput = [];
        trialIdx = {spkData.(thisSesActiveUnits{un,1}{1,1}).TrialIdx}';
        unitActiveTrials = {spkData.(thisSesActiveUnits{un,1}{1,1}).UnitActiveTrials}';
        alignedSpikeTimes = {spkData.(thisSesActiveUnits{un,1}{1,1}).AlignedSpikeTimes}';
        spikeCount = aa;
        clear aa
        for ii = 1:size(uniqEVandVar,1)
            iidx = find(ismember(bb(:,2:3),uniqEVandVar(ii,:),'rows'));
            tempOutput = [ tempOutput; struct('TrialIdx', [trialIdx{iidx(1)}; trialIdx{iidx(2)}], ...
                                              'Value', uniqEVandVar(ii,:), ...
                                              'UnitActiveTrials',  [unitActiveTrials{iidx(1)}; unitActiveTrials{iidx(2)}], ...
                                              'SpikeCount',  [spikeCount{iidx(1)}; spikeCount{iidx(2)}]) ] ;
            clear iidx                                              
        end
        clear ii bb fieldNames trialIdx unitActiveTrials alignedSpikeTimes alignedSpikeTimes spikeCount uniqEVandVar
        spkData.(thisSesActiveUnits{un,1}{1,1}) = tempOutput;
        clear tempOutput

        switch [num2str(thisActiveChanArrayIdx) ' ' num2str(~isempty(spkData.(thisSesActiveUnits{un,1}{1,1})))]
          case '1 1'
            popUnitIdxArr1 = popUnitIdxArr1 + 1;
            sz = cellfun(@(x) length(x), {spkData.(thisSesActiveUnits{un,1}{1,1}).SpikeCount}');
            thisUnitParaVal = cellfun( @(x,y) repmat(x, y, 1), {spkData.(thisSesActiveUnits{un,1}{1,1}).Value}', mat2cell(sz, ones(length(sz),1) ,1) , 'UniformOutput' , 0 ); 
            thisUnitParaVal = cell2mat(thisUnitParaVal);
            clear sz
            
            % collect firing rate of the neuron, trial-by-trial
            thisUnitSpkCount = cell2mat({spkData.(thisSesActiveUnits{un,1}{1,1}).SpikeCount}'); % each row will be the spike count for one unit
            % zscore firing rates of this neuron, accross trials, take care of the NaN values
            aa = zscore(thisUnitSpkCount(~isnan(thisUnitSpkCount)));
            thisUnitSpkCount(~isnan(thisUnitSpkCount)) = aa;
            clear aa
            
            % add this unit's trial-by-trial firing rate, parameter values for each trila and unit index to the population
            allUnitTrlsSpkCountArr1 = [allUnitTrlsSpkCountArr1; thisUnitSpkCount];
            allUnitTrlsParaValArr1  = [allUnitTrlsParaValArr1;  thisUnitParaVal ];
            allUnitTrlsUnitIdxArr1  = [allUnitTrlsUnitIdxArr1;  repmat(popUnitIdxArr1, size(thisUnitParaVal,1),1)];
            allUnitGroupsArr1       = [allUnitGroupsArr1; { {spkData.(thisSesActiveUnits{un,1}{1,1}).SpikeCount}' }];
            
          case '2 1'
            popUnitIdxArr2 = popUnitIdxArr2 + 1;
            sz = cellfun(@(x) length(x), {spkData.(thisSesActiveUnits{un,1}{1,1}).SpikeCount}');
            thisUnitParaVal = cellfun( @(x,y) repmat(x, y, 1), {spkData.(thisSesActiveUnits{un,1}{1,1}).Value}', mat2cell(sz, ones(length(sz),1) ,1) , 'UniformOutput' , 0 ); 
            thisUnitParaVal = cell2mat(thisUnitParaVal);
            clear sz
            
            % collect firing rate of the neuron, trial-by-trial
            thisUnitSpkCount = cell2mat({spkData.(thisSesActiveUnits{un,1}{1,1}).SpikeCount}'); % each row will be the spike count for one unit
            % zscore firing rates of this neuron, accross trials, take care of the NaN values
            aa = zscore(thisUnitSpkCount(~isnan(thisUnitSpkCount)));
            thisUnitSpkCount(~isnan(thisUnitSpkCount)) = aa;
            clear aa
            
            % add this unit's trial-by-trial firing rate to the population
            allUnitTrlsSpkCountArr2 = [allUnitTrlsSpkCountArr2; thisUnitSpkCount];
            allUnitTrlsParaValArr2  = [allUnitTrlsParaValArr2;  thisUnitParaVal ];
            allUnitTrlsUnitIdxArr2  = [allUnitTrlsUnitIdxArr2;  repmat(popUnitIdxArr2, size(thisUnitParaVal,1),1)];
            allUnitGroupsArr2       = [allUnitGroupsArr2; { {spkData.(thisSesActiveUnits{un,1}{1,1}).SpikeCount}' }];
            
            clear thisUnitParaVal thisUnitSpkCount
        end
    end
    clear thisSesActiveUnits un spkData*
    ses = ses + 1;
end
clear ses popUnitIdxArr1 popUnitIdxArr2

% remove all NaN trials from the population data
nanIdx1 = isnan(allUnitTrlsSpkCountArr1);
allUnitTrlsSpkCountArr1(nanIdx1) = [];
allUnitTrlsParaValArr1(nanIdx1,:) = [];    
allUnitTrlsUnitIdxArr1(nanIdx1) = [];        
    
nanIdx2 = isnan(allUnitTrlsSpkCountArr2);
allUnitTrlsSpkCountArr2(nanIdx2) = [];
allUnitTrlsParaValArr2(nanIdx2,:) = [];    
allUnitTrlsUnitIdxArr2(nanIdx2) = [];        

clear nanIdx*
%%  Get the discriminability index across conditions for individual units
unitDiscrMatrix1 = cell(length(unique(allUnitTrlsUnitIdxArr1)),1); % length(unitDiscrMatrixArray1) = number of units on array1 
unitDiscrMatrix2 = cell(length(unique(allUnitTrlsUnitIdxArr2)),1); % length(unitDiscrMatrixArray2) = number of units on array2 

% initialize the discriminablity matrix
unitDiscrMatrix1 = repmat({nan(size(allUnitGroupsArr1{1,1},1))}, size(allUnitGroupsArr1)); 
unitDiscrMatrix2 = repmat({nan(size(allUnitGroupsArr2{1,1},1))}, size(allUnitGroupsArr2)); 

    switch DiscriminabilityType 
        case 'Cohen''s d'
            unitDiscrMatrix1 = cellfun(@(u)  cell2mat( cellfun(@(x)  cellfun(@(x1,y)  (nanmean(x1)-nanmean(y))/sqrt( 0.5*((nanstd(x1)^2)+(nanstd(y)^2)) ) , repmat({x},size(u,1),1), u) ,  u, 'UniformOutput', 0)'), ...
                                       allUnitGroupsArr1, 'UniformOutput', 0);

            unitDiscrMatrix2 = cellfun(@(u)  cell2mat( cellfun(@(x)  cellfun(@(x1,y)  (nanmean(x1)-nanmean(y))/sqrt( 0.5*((nanstd(x1)^2)+(nanstd(y)^2)) ) , repmat({x},size(u,1),1), u) ,  u, 'UniformOutput', 0)'), ...
                                       allUnitGroupsArr2, 'UniformOutput', 0);
    end
% % % % unitDiscrMatrix1 = cellfun(@(u)  cell2mat( cellfun(@(x)  cellfun(@(x1,y)  (nanmean(x1)-nanmean(y))/sqrt((nanvar(x1)/length(x1))+(nanvar(y)/length(y))) , repmat({x},size(u,1),1), u) ,  u, 'UniformOutput', 0)'), ...
% % % %                            allUnitGroupsArr1, 'UniformOutput', 0);
% % % % 
% % % % unitDiscrMatrix2 = cellfun(@(u)  cell2mat( cellfun(@(x)  cellfun(@(x1,y)  (nanmean(x1)-nanmean(y))/sqrt((nanvar(x1)/length(x1))+(nanvar(y)/length(y))) , repmat({x},size(u,1),1), u) ,  u, 'UniformOutput', 0)'), ...
% % % %                            allUnitGroupsArr2, 'UniformOutput', 0);

%% plot lick vs units discriminability
fg1 = figure('Name', [Monkey ', Lick vs individual cells, ' DiscriminabilityType])
numRows = 2; % number of areas
numColumns = 2; % number of time intervals for lick signal
plotCounter = 0;
for lickTimIntrvl = 1:length(lickDiscrMatrix)
    plotCounter = plotCounter+1;
    % ---- dlPFC ----
    subplot(numRows, numColumns, plotCounter), hold on, box on
    titl = [ {['dlPFC']} ,...
             {['lick ' alignEvent{lickTimIntrvl} '[' num2str(eventTimeInterval{lickTimIntrvl}(1)) ' , ' num2str(eventTimeInterval{lickTimIntrvl}(2)) ']' ]} ,...
             {['spikes: ' spikeAlignedEvent '[' num2str(spikeTimeInterval(1)) ' , ' num2str(spikeTimeInterval(2)) ']' ]} ];
    title(titl)
    xlabel('lick discriminability index')
    ylabel('units discriminability index')
       
    for cd1 = 1 : size(lickDiscrMatrix{lickTimIntrvl},1)-1
        for cd2 = cd1+1 : size(lickDiscrMatrix{lickTimIntrvl},1)
            cellfun(@(x) scatter(abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) , abs(x(cd1,cd2)), 'Marker', 'o', 'MarkerEdgeColor', 'b') , ...
                    unitDiscrMatrix1, 'UniformOutput', 0)
        end
        clear cd2
    end
    clear cd1 titl
    plot([min([xlim ylim]) max([xlim ylim])],[min([xlim ylim]) max([xlim ylim])],'--', 'Color', 'k')
    set(gca, 'XLim', [min([xlim ylim]) max([xlim ylim])], 'YLim', [min([xlim ylim]) max([xlim ylim])])
    axis square 


    % ---- 7a ----
    subplot(numRows, numColumns, plotCounter+numColumns), hold on, box on
    titl = [ {['7a']} ,...
             {['lick ' alignEvent{lickTimIntrvl} '[' num2str(eventTimeInterval{lickTimIntrvl}(1)) ' , ' num2str(eventTimeInterval{lickTimIntrvl}(2)) ']' ]} ,...
             {['spikes: ' spikeAlignedEvent '[' num2str(spikeTimeInterval(1)) ' , ' num2str(spikeTimeInterval(2)) ']' ]} ];
    title(titl)
    xlabel('lick discriminability index')
    ylabel('units discriminability index')
       
    for cd1 = 1 : size(lickDiscrMatrix{lickTimIntrvl},1)-1
        for cd2 = cd1+1 : size(lickDiscrMatrix{lickTimIntrvl},1)
            cellfun(@(x) scatter(abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) , abs(x(cd1,cd2)), 'Marker', 'o', 'MarkerEdgeColor', 'r') , ...
                    unitDiscrMatrix2, 'UniformOutput', 0)
        end
        clear cd2
    end
    clear cd1 titl
    plot([min([xlim ylim]) max([xlim ylim])],[min([xlim ylim]) max([xlim ylim])],'--', 'Color', 'k')
    set(gca, 'XLim', [min([xlim ylim]) max([xlim ylim])], 'YLim', [min([xlim ylim]) max([xlim ylim])])
    axis square 
    
end

%% plot distribution of regression line slope for individual unit vs lick(all sessions together!) discriminability 
% this is the indices of the lower triangular part of the discriminability matrix
indx1 = find(tril(ones(size(unitDiscrMatrix1{1})),-1));
indx2 = find(tril(ones(size(unitDiscrMatrix2{1})),-1));

unitsDisVector1 = cellfun(@(x) abs(x(indx1)), unitDiscrMatrix1, 'UniformOutput', 0);
unitsDisVector2 = cellfun(@(x) abs(x(indx2)), unitDiscrMatrix2, 'UniformOutput', 0);

fg2 = figure('Name', [Monkey ', distribution of regression line between units and lick discriminability, ' DiscriminabilityType], 'Position', [2308 460 1262 420]);
plotCounter = 0;

for tim = 1:length(lickDiscrMatrix)
    [regSlope1{tim},regConfInterval1{tim},regResidual1{tim}] = cellfun(@(y) regress( y, abs(lickDiscrMatrix{tim}(indx1))) , unitsDisVector1 ,'UniformOutput', 0);
    [regSlope2{tim},regConfInterval2{tim},regResidual2{tim}] = cellfun(@(y) regress( y, abs(lickDiscrMatrix{tim}(indx2))) , unitsDisVector2 ,'UniformOutput', 0);

    plotCounter = plotCounter+1;
    subplot(length(lickDiscrMatrix),2,plotCounter), hold on, box on
    xlabel('slope of regression line of spike vs lick discriminability')
    ylabel('proportion of units')
    titl = [ {[Monkey ', dlPFC']} ,...
             {['lick ' alignEvent{1} '[' num2str(eventTimeInterval{tim}(1)) ' , ' num2str(eventTimeInterval{tim}(2)) ']' ]} ,...
             {['spikes: ' spikeAlignedEvent '[' num2str(spikeTimeInterval(1)) ' , ' num2str(spikeTimeInterval(2)) ']' ]} ];
    title(titl)
    clear titl 
    binSize = 0.05;
    binLeftEdge = min(0,min(cell2mat(regSlope1{tim}))-(2*binSize)) : binSize : max(cell2mat(regSlope1{tim}))+(2*binSize);
    N = histcounts(cell2mat(regSlope1{tim, 1}),binLeftEdge);
    bar(binLeftEdge(1:end-1)+(0.5*binSize), N/sum(N), 'FaceColor', [0 0.44 0.74]);
    clear N

    plotCounter = plotCounter+1;
    subplot(length(lickDiscrMatrix),2,plotCounter), hold on, box on
    xlabel('slope of regression line of spike vs lick discriminability')
    ylabel('proportion of units')
    titl = [ {[Monkey ', 7a']} ,...
             {['lick ' alignEvent{1} '[' num2str(eventTimeInterval{tim}(1)) ' , ' num2str(eventTimeInterval{tim}(2)) ']' ]} ,...
             {['spikes: ' spikeAlignedEvent '[' num2str(spikeTimeInterval(1)) ' , ' num2str(spikeTimeInterval(2)) ']' ]} ];
    title(titl)
    clear titl       
    N = histcounts(cell2mat(regSlope2{tim, 1}),binLeftEdge);
    bar(binLeftEdge(1:end-1)+(0.5*binSize), N/sum(N), 'FaceColor', [0 0.44 0.74]);
    clear N binLeftEdge binSize

    saveas(fg2, [save_dir Monkey '_HistRegresSlope_SpikeLick_' num2str(tim)], 'fig')
    saveas(fg2, [save_dir Monkey '_HistRegresSlope_SpikeLick_' num2str(tim)], 'bmp')
    close(fg2)
    clear fg2
end
clear tim indx* tim

%%
aa1 = reshape( cell2mat(unitDiscrMatrix1'), [size(unitDiscrMatrix1{1})  size(unitDiscrMatrix1,1)]); % size(aa) = [#conditions #conditions #unitsOnThisArray]
fg3 = figure('Name', [Monkey ', average units discriminability, ' DiscriminabilityType])
subplot(1,2,1), hold on
title('dlPFC')
for i =1:size(aa1,1)-1
    for j=i+1:size(aa1,2)
        plot3(i, j, nanmean(abs(aa1(i,j,:))), 'o')
        plot3([i i], [j j], [nanmean(abs(aa1(i,j,:))) - nanstd(abs(aa1(i,j,:))) , nanmean(abs(aa1(i,j,:))) + nanstd(abs(aa1(i,j,:)))])
    end
end
clear i j
surf(nanmean(abs(aa1),3))

aa2 = reshape( cell2mat(unitDiscrMatrix2'), [size(unitDiscrMatrix2{1})  size(unitDiscrMatrix2,1)]); % size(aa) = [#conditions #conditions #unitsOnThisArray]
subplot(1,2,2), hold on
title('7a')
for i =1:size(aa2,1)-1
    for j=i+1:size(aa2,2)
        plot3(i, j, nanmean(abs(aa2(i,j,:))), 'o')
        plot3([i i], [j j], [nanmean(abs(aa2(i,j,:))) - nanstd(abs(aa2(i,j,:))) , nanmean(abs(aa2(i,j,:))) + nanstd(abs(aa2(i,j,:)))])
    end
end
clear i j
surf(nanmean(abs(aa2),3))

%%
colorCode = [0 1 1
             0 0 1
             0 0 0
             ...
             1 0.6 0.6    
             0.99 0.2 0.2
             0.6 0 0
             ...
             0.6 1 0.6
             0 1 0
             0.1 0.5 0.1];

matrix_unitDiscrMatrix1 = reshape( cell2mat(unitDiscrMatrix1'), [size(unitDiscrMatrix1{1})  size(unitDiscrMatrix1,1)]); % size(.) = [#conditions #conditions #unitsOnThisArray]
matrix_unitDiscrMatrix2 = reshape( cell2mat(unitDiscrMatrix2'), [size(unitDiscrMatrix2{1})  size(unitDiscrMatrix2,1)]); % size(.) = [#conditions #conditions #unitsOnThisArray]


% get the parameter values for each condition
grpParaVal = {allSesGroupedData{1,1}.Value}';
% grpParaValString = cellfun(@(x) ['[' num2str(x(1)) ',' num2str(x(2)) ']'], grpParaVal, 'UniformOutput', 0);

% fg2 = figure('Name', [Monkey ', 2 areas'])
fg4_dlpfc = figure('Name', [Monkey ', dlPFC, ' DiscriminabilityType])
fg4_7a = figure('Name', [Monkey ', 7a, ' DiscriminabilityType])

numRows = 1;
numColumns = 2;
plotCounter = 0;
for lickTimIntrvl = 1:length(lickDiscrMatrix)
       plotCounter = plotCounter+1; 
       
       % scatter plot dlPFC
       figure(fg4_dlpfc)
       subplot(numRows, numColumns, plotCounter), hold on, box on
       titl = [ {['dlPFC']} ,...
                {['lick ' alignEvent{lickTimIntrvl} '[' num2str(eventTimeInterval{lickTimIntrvl}(1)) ' , ' num2str(eventTimeInterval{lickTimIntrvl}(2)) ']' ]} ,...
                {['spikes: ' spikeAlignedEvent '[' num2str(spikeTimeInterval(1)) ' , ' num2str(spikeTimeInterval(2)) ']' ]} ];
       title(titl)
       xlabel('lick discriminability index')
       ylabel('average units discriminability index')
       xVal = [];
       yVal = [];
       for cd1 = 1 : size(lickDiscrMatrix{lickTimIntrvl},1)-1
           for cd2 = cd1+1 : size(lickDiscrMatrix{lickTimIntrvl},1)
               % plot dlPFC
               plot(abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) , nanmean(abs(matrix_unitDiscrMatrix1(cd1,cd2,:))), 'Marker', 'o', ...
                   'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerEdgeColor', colorCode(cd1,:), 'MarkerFaceColor', colorCode(cd2,:))
               xVal = [xVal; abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2))];
               yVal = [yVal; nanmean(abs(matrix_unitDiscrMatrix1(cd1,cd2,:)))];
           end
           clear cd2
       end
       clear cd1 titl       
       [b,bint] = regress(yVal,[ones(size(xVal)) xVal]);
       plot(xlim, b(1)+(b(2)*xlim), 'k');
       text(mean(xlim), mean(ylim), [{['Slope ' num2str(b(2))]},{['conf. interval [' num2str(bint(2,1)) ',' num2str(bint(2,2)) ']']}])
       clear b bint xVal yVal
       
       % scatter plot 7a
       figure(fg4_7a)
       subplot(numRows, numColumns, plotCounter), hold on, box on
       titl = [ {['7a']} ,...
                {['lick ' alignEvent{lickTimIntrvl} '[' num2str(eventTimeInterval{lickTimIntrvl}(1)) ' , ' num2str(eventTimeInterval{lickTimIntrvl}(2)) ']' ]} ,...
                {['spikes: ' spikeAlignedEvent '[' num2str(spikeTimeInterval(1)) ' , ' num2str(spikeTimeInterval(2)) ']' ]} ];
       title(titl)
       xlabel('lick discriminability index')
       ylabel('average units discriminability index')
       xVal = [];
       yVal = [];
       for cd1 = 1 : size(lickDiscrMatrix{lickTimIntrvl},1)-1
           for cd2 = cd1+1 : size(lickDiscrMatrix{lickTimIntrvl},1)
               % plot 7a
               plot(abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2)) , nanmean(abs(matrix_unitDiscrMatrix2(cd1,cd2,:))), 'Marker', 'o', ...
                   'LineWidth', 2.5, 'MarkerSize', 10, 'MarkerEdgeColor', colorCode(cd1,:), 'MarkerFaceColor', colorCode(cd2,:))
               xVal = [xVal; abs(lickDiscrMatrix{lickTimIntrvl}(cd1,cd2))];
               yVal = [yVal; nanmean(abs(matrix_unitDiscrMatrix2(cd1,cd2,:)))];
           end
           clear cd2
       end
       clear cd1 titl        
       [b,bint] = regress(yVal,[ones(size(xVal)) xVal]);
       plot(xlim, b(1)+(b(2)*xlim), 'k');
       text(mean(xlim), mean(ylim), [{['Slope ' num2str(b(2))]},{['conf. interval [' num2str(bint(2,1)) ',' num2str(bint(2,2)) ']']}])
       clear b bint xVal yVal
       
       
       % line plot
       if 0
        plot(1:size(pairedCdIdx,1) , [sign(diff(pairedCdDiscrimIdxLick)) 0]  , 'Color', 'r')
        plot(1:size(pairedCdIdx,1) , [sign(diff(pairedCdDiscrimIdxLFPArr1)) 0], 'Color', 'c')
        plot(1:size(pairedCdIdx,1) , [sign(diff(pairedCdDiscrimIdxLFPArr2)) 0], 'Color', 'b')
       
        set(gca, 'XTick', 1:size(pairedCdIdx,1), 'XTickLabel', xTickLabel)
        xtickangle(gca,90)
        legend('Lick','dlPFC','7a')
       end
       
end

saveas(fg4_dlpfc, [save_dir Monkey '_Scatter_MeanUnitsDiscVsLick_dlPFC'], 'fig')
saveas(fg4_dlpfc, [save_dir Monkey '_Scatter_MeanUnitsDiscVsLick_dlPFC'], 'bmp')

saveas(fg4_7a, [save_dir Monkey '_Scatter_MeanUnitsDiscVsLick_7a'], 'fig')
saveas(fg4_7a, [save_dir Monkey '_Scatter_MeanUnitsDiscVsLick_7a'], 'bmp')

close(fg4_dlpfc), clear fg4_dlpfc
close(fg4_7a), clear fg4_7a
