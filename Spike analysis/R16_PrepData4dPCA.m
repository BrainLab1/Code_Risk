% make input data for dPCA

% 08.05.2018: this is the same scripts as R13_PrepData4dPCA but with the new channel-electrode map 
% created 10.03.2018
% -------------------------------------------------------------------------

clear, clc
dbstop if error 

%%
spike_dir   = 'X:\data\Add PSTH to Aligned Spk\'; % 

Monkey = 'MacDuff'; % select among 'Mojo', 'MacDuff' and 'Data'
                            % if you would like to select both monkeys, use
                            % 'Data' as the monkey name
alignedEvent = '(cue)';  % select among 'cue', 'reward' and    % ***** also add the option to plot the overall psth, not grouped by condition
                                        % 'Data'. If you would like to
                                        % select both events, use 'Data'
                                        % as the event name 

Area = '7a';
ArrayIdx = find(strcmp({'dlPFC','7a'},Area));

PoolUnitsFromAllSessions = 1;

glmPara = 'ev'; % 'ev'; 'var'; 'int'
timeIntervalOfInterest = [0, 0.8];      
sigTimeInterval = [0.2 0.4];
minRepetitionPerCondition = 5;
groupingCriteriaForDataPath = 'CuePos & CurrEV & CurrVar';
groupingCriteria = 'expecter_reward & RewardVariance'; %'Variance'; %'expecter_reward & RewardVariance';
numNoise= -10;
% % % % SigUnitPath = ['X:\Report 13\Spike Counts\' alignedEvent(2:end-1) '[' num2str(sigTimeInterval(1)) ',' num2str(sigTimeInterval(2)) ']_' groupingCriteriaForDataPath '_minRepetition' num2str(minRepetitionPerCondition) '\'];

% we need these parameters to read out the right file for list of the active units. these values are included in the file names.
StabilityCheckTimeInterval = [-0.5 1];
AveSpkCountThreshold = 5;
ActiveUnitsPath = ['X:\data\' Monkey '_UnitStability_' alignedEvent '[' regexprep(num2str(StabilityCheckTimeInterval(1)),'0.','0,') ' ' regexprep(num2str(StabilityCheckTimeInterval(2)),'0.','0,') ']_AveSpkCntThreshold' num2str(AveSpkCountThreshold)];

RemoveSigUnitsFromAnalysis = 0;
CorrectForBaseLine = 1;

%%  load significant units list if needed
% % % % if RemoveSigUnitsFromAnalysis
% % % %     load([SigUnitPath 'Results\Active Units\' Monkey '_SigUnitIDs_' Area '_cue'])
% % % %     % this will load following variables 
% % % %     %     'allUnitIDsArray1','evSigUnitIdx','varSigUnitIdx','intSigUnitIdx','pGLMUnitArr','bGLMUnitArr',
% % % %     %     'timeInterval','alignedEvent', 'allUnitTrlsSpkCountArr', 'allUnitTrlsUnitIdxArr'
% % % % 
% % % %     % read out session names for which there has been signficant units
% % % %     eval(['sigSessions = unique(allUnitIDsArray' num2str(ArrayIdx) '(' glmPara 'SigUnitIdx, 1));'])
% % % % 
% % % %     % read out significant session/unit names
% % % %     eval(['sigUnitIDsArray = allUnitIDsArray' num2str(ArrayIdx) '(' glmPara 'SigUnitIdx, :);'])
% % % % % % % % 
% % % % % % % %     % read out GLM betta values for significant units
% % % % % % % %     eval(['bSigUnits = bGLMUnitArr' num2str(ArrayIdx) '(' glmPara 'SigUnitIdx,:);'])   
% % % % end
%% load active units from both arrays; re-assigne units to arrays using new electrod map; get the active units list for this array
% load list of all active units
load(ActiveUnitsPath, 'ActiveSesUnitsArr1' , 'ActiveSesUnitsArr2' );
% get the new channel-electrode assignments on the array
eval(['[newArray1, newArray2] = new_electrodepinout' Monkey '(''chan'');']) 

ActiveSesUnitsArr1(:,2) = cellfun(@(x) x{1,1}, ActiveSesUnitsArr1(:,2), 'UniformOutput', 0);
ActiveSesUnitsArr2(:,2) = cellfun(@(x) x{1,1}, ActiveSesUnitsArr2(:,2), 'UniformOutput', 0);

charIdx1 = strfind(ActiveSesUnitsArr1(:,2),'ch_');
charIdx2 = strfind(ActiveSesUnitsArr2(:,2),'ch_');

chanNum1 = cellfun(@(x,y) str2num(x(y+3:end)), ActiveSesUnitsArr1(:,2), charIdx1);
chanNum2 = cellfun(@(x,y) str2num(x(y+3:end)), ActiveSesUnitsArr2(:,2), charIdx2);

newActiveSesUnitsArr1 = [ActiveSesUnitsArr1(find(ismember(chanNum1,newArray1)),:) ; ActiveSesUnitsArr2(find(ismember(chanNum2,newArray1)),:)];
newActiveSesUnitsArr2 = [ActiveSesUnitsArr1(find(ismember(chanNum1,newArray2)),:) ; ActiveSesUnitsArr2(find(ismember(chanNum2,newArray2)),:)];

eval(['ActiveSesUnitsArr = newActiveSesUnitsArr' num2str(ArrayIdx) ';'])
eval(['sesListPSTH = unique(ActiveSesUnitsArr(:,1));'])

clear ActiveSesUnitsArr1 ActiveSesUnitsArr2 newActiveSesUnitsArr1 newActiveSesUnitsArr2

% % % % eval(['load(ActiveUnitsPath, ''ActiveSesUnitsArr' num2str(ArrayIdx) ''' );'])
% % % % eval(['ActiveSesUnitsArr = ActiveSesUnitsArr' num2str(ArrayIdx) ';'])
% % % % eval(['clear ActiveSesUnitsArr' num2str(ArrayIdx) ])
% % % % eval(['sesListPSTH = unique(ActiveSesUnitsArr(:,1));'])

%% 
idx = 1;
% read out the session files one by one
while idx <= length(sesListPSTH)
    display(['ses ' num2str(idx)])
    % load the spikeTrials for current sessions
    load ([spike_dir sesListPSTH{idx}]);

    % read out the PSTH time axes
    T = spikeTrials.PSTH.psthTime;   % ********  add the option to plot the psth in desired time interval

    % get the indices for time intreval of interest
    [~, t1Idx] = min( abs(timeIntervalOfInterest(1)-T) );
    [~, t2Idx] = min( abs(timeIntervalOfInterest(2)-T) );
    
    % find active units on array 1 for this session
    thisSesActiveUnits = ActiveSesUnitsArr( find( strcmp(ActiveSesUnitsArr(:,1),sesListPSTH{idx}) ) , 2);
% % % %     thisSesActiveUnits = cellfun(@(x)   x{1,1}, thisSesActiveUnits, 'UniformOutput', 0);
    
    % remove significant units from the list of active units for array 1
    if RemoveSigUnitsFromAnalysis
        % find significant units in this session if asked
        thisSesSigUnitIdx = find(strcmp(sesListPSTH{idx},sigUnitIDsArray(:,1)));
        idxToRemove = cellfun(@(x) find(strcmp(thisSesActiveUnits, x)) , sigUnitIDsArray(thisSesSigUnitIdx,2) , 'UniformOutput' , 0 );
        thisSesActiveUnits(cell2mat(idxToRemove),:)=[];
    end
    
    % get the PSTH for active units on the array
    for un = 1:length(thisSesActiveUnits)       
        % zscore the unit psth within the session
        aa = spikeTrials.PSTH.(thisSesActiveUnits{un}).psth;
        bb = reshape(aa',1,[]);
        cc = zscore(bb);
        zUnitPSTH = reshape(cc , size(aa'))';
        clear aa bb cc  

        [output] = FindRepetitionPerConditionForUnit (groupingCriteria, spikeTrials.cfg.event, spikeTrials.PSTH.(thisSesActiveUnits{un}));
        % initialize variables 
        if ~exist('allUnitsPSTH_Arr','var')
            allUnitsPSTH_Arr  = cell(length(output),length(sesListPSTH)); % cell{grp, session} = matrix: each row is the averaged PSTH for each unit; 
            allNoiseUnitPSTH_Arr = cell(length(output),length(sesListPSTH), numNoise); % cell{grp, session, noiseIdx} = matrix: each row is the averaged PSTH for each unit; 
            groupingInfo.GroupingParaVal = {output.Value}';
        end
        
        for grp = 1:length(output)
            % find for which trials in this group, the unit was active
            unitValidTrlsInGrp = output(grp).TrialIdx(output(grp).UnitActiveTrials);
            % find which rows of the psth matrix are relevant for this group
            unitValidPSTHinGrp = ismember(spikeTrials.PSTH.(thisSesActiveUnits{un}).active_eventStrctTrls , unitValidTrlsInGrp);
            % get mean psth for this group 
            allUnitsPSTH_Arr{grp, idx} = [allUnitsPSTH_Arr{grp, idx}; mean( zUnitPSTH(unitValidPSTHinGrp, t1Idx:t2Idx) , 1 )];
            clear unitValidPSTHinGrp unitValidTrlsInGrp
        end
        clear grp output 
        
        % randomise trials in event structure to generate noise
        ns = 0;
        while ns < numNoise
            ns = ns + 1;
            [un ns]
            unitActiveTrls_event = spikeTrials.PSTH.(thisSesActiveUnits{un}).active_eventStrctTrls;
            unitEventStrct = spikeTrials.cfg.event(unitActiveTrls_event);
            noiseTrlNum = randperm(length(unitActiveTrls_event))';
            noiseEventStrct = spikeTrials.cfg.event;
            noiseEventStrct(unitActiveTrls_event) = unitEventStrct(noiseTrlNum);
            clear unitEventStrct noiseTrlNum unitActiveTrls_event
            [outputNoise{ns}] = FindRepetitionPerConditionForUnit (groupingCriteria, noiseEventStrct, spikeTrials.PSTH.(thisSesActiveUnits{un}));
            for grp = 1:length(outputNoise{ns})
                % find for which trials in this group, the unit was active
                unitValidTrlsInGrp = outputNoise{ns}(grp).TrialIdx(outputNoise{ns}(grp).UnitActiveTrials);
                % find which rows of the psth matrix are relevant for this group
                unitValidPSTHinGrp = ismember(spikeTrials.PSTH.(thisSesActiveUnits{un}).active_eventStrctTrls , unitValidTrlsInGrp);
                % get mean psth for this group 
                allNoiseUnitPSTH_Arr{grp, idx, ns} = [allNoiseUnitPSTH_Arr{grp, idx, ns}; mean( zUnitPSTH(unitValidPSTHinGrp, t1Idx:t2Idx) , 1 )];
                clear unitValidPSTHinGrp unitValidTrlsInGrp
            end
            clear grp output 
        end
        clear ns zUnitPSTH
    end
    
    clear un thisSesActiveUnits unitNamesArr2 spikeTrials thisSesSigUnitIdx
    idx = idx + 1; % go to the next session
end % end of sessions
clear idx 

newT = T(t1Idx:t2Idx);

% baseline correction -----
if CorrectForBaseLine
    [~ , t2BaseIdx] = min(abs(newT));
    % % % % t2BaseIdx = 1;
    allUnitsPSTH_Arr_Base = cellfun(@(x) x - mean(x(:,1:t2BaseIdx),2)  , allUnitsPSTH_Arr, 'UniformOutput', 0);
    clear t2BaseIdx
else 
    allUnitsPSTH_Arr_Base = allUnitsPSTH_Arr;
end

%% pool units from all sessions

if PoolUnitsFromAllSessions
    aa = allUnitsPSTH_Arr_Base;
    for i=1:size(aa,1)
        bb{i,1}=cell2mat(aa(i,:)');
    end
    allUnitsPSTH_Arr_Base = bb;
    clear bb i
end
