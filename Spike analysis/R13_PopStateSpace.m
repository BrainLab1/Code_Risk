
% This script makes the high dimensional neural space and plots trajectory of neural activity for different task 
% conditions in the PCA space; also plots PCs' PSTHs 

% created 07.03.2018
% -------------------------------------------------------------------------

clear, clc
dbstop if error 

%%
spike_dir   = 'X:\data\Add PSTH to Aligned Spk\'; % 

% % % % colorCode = [0.6 0.6 1; 0.2 0.3 1; 0 0 0;...
% % % %               0.2 0.6 1; 0.4 0.3 1; 0.6 0 0;...
% % % %               0.6 1 0.6; 0.9 0.2 1; 0.6 0 0.5];   % ***** use specific color for each task condition; 36 colors

Monkey = 'MacDuff'; % select among 'Mojo', 'MacDuff' and 'Data'
                            % if you would like to select both monkeys, use
                            % 'Data' as the monkey name
alignedEvent = '(cue)';  % select among 'cue', 'reward' and    % ***** also add the option to plot the overall psth, not grouped by condition
                                        % 'Data'. If you would like to
                                        % select both events, use 'Data'
                                        % as the event name 

Area = 'dlPFC';
ArrayIdx = find(strcmp({'dlPFC','7a'},Area));

PoolUnitsFromAllSessions = 1;

glmPara = 'ev'; % 'ev'; 'var'; 'int'
timeIntervalOfInterest = [0, 0.8];      
sigTimeInterval = [0 0.2];
minRepetitionPerCondition = 5;
groupingCriteriaForDataPath = 'CuePos & CurrEV & CurrVar';
groupingCriteria = 'expecter_reward & RewardVariance'; %'Variance'; %'expecter_reward & RewardVariance';
numNoise= -10;
SigUnitPath = ['X:\Report 13\Spike Counts\' alignedEvent(2:end-1) '[' num2str(sigTimeInterval(1)) ',' num2str(sigTimeInterval(2)) ']_' groupingCriteriaForDataPath '_minRepetition' num2str(minRepetitionPerCondition) '\'];

% we need these parameters to read out the right file for list of the active units. these values are included in the file names.
StabilityCheckTimeInterval = [-0.5 1];
AveSpkCountThreshold = 5;
ActiveUnitsPath = ['X:\data\' Monkey '_UnitStability_' alignedEvent '[' regexprep(num2str(StabilityCheckTimeInterval(1)),'0.','0,') ' ' regexprep(num2str(StabilityCheckTimeInterval(2)),'0.','0,') ']_AveSpkCntThreshold' num2str(AveSpkCountThreshold)];

RemoveSigUnitsFromAnalysis = 0;

%%  load significant units list if needed
if RemoveSigUnitsFromAnalysis
    load([SigUnitPath 'Results\Active Units\' Monkey '_SigUnitIDs_' Area '_cue'])
    % this will load following variables 
    %     'allUnitIDsArray1','evSigUnitIdx','varSigUnitIdx','intSigUnitIdx','pGLMUnitArr','bGLMUnitArr',
    %     'timeInterval','alignedEvent', 'allUnitTrlsSpkCountArr', 'allUnitTrlsUnitIdxArr'

    % read out session names for which there has been signficant units
    eval(['sigSessions = unique(allUnitIDsArray' num2str(ArrayIdx) '(' glmPara 'SigUnitIdx, 1));'])

    % read out significant session/unit names
    eval(['sigUnitIDsArray = allUnitIDsArray' num2str(ArrayIdx) '(' glmPara 'SigUnitIdx, :);'])
% % % % 
% % % %     % read out GLM betta values for significant units
% % % %     eval(['bSigUnits = bGLMUnitArr' num2str(ArrayIdx) '(' glmPara 'SigUnitIdx,:);'])   
end
%% get the active units list for this array
eval(['load(ActiveUnitsPath, ''ActiveSesUnitsArr' num2str(ArrayIdx) ''' );'])
eval(['ActiveSesUnitsArr = ActiveSesUnitsArr' num2str(ArrayIdx) ';'])
eval(['clear ActiveSesUnitsArr' num2str(ArrayIdx) ])
eval(['sesListPSTH = unique(ActiveSesUnitsArr(:,1));'])

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
    thisSesActiveUnits = cellfun(@(x)   x{1,1}, thisSesActiveUnits, 'UniformOutput', 0);
    
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
[~ , t2BaseIdx] = min(abs(newT));
% % % % t2BaseIdx = 1;
allUnitsPSTH_Arr_Base = cellfun(@(x) x - mean(x(:,1:t2BaseIdx),2)  , allUnitsPSTH_Arr, 'UniformOutput', 0);
clear t2BaseIdx


% % % % save('X:\Report 12\allUnitsPSTH_Arr_Base','allUnitsPSTH_Arr_Base')
% % % % save('X:\Report 12\allUnitsPSTH_Arr2_Base','allUnitsPSTH_Arr2_Base')
% % % % save('X:\Report 12\allNoiseUnitPSTH_Arr','allNoiseUnitPSTH_Arr')
% % % % save('X:\Report 12\allNoiseUnitPSTH_Arr2','allNoiseUnitPSTH_Arr2')
% % % % save('X:\Report 12\Workspace')
%% test trajectories' distance in the high dimentional neural space
if 0
ses = 11;
pairConditions = nchoosek([1:length(groupingInfo.GroupingParaVal)],2); 

for i = 1:size(pairConditions,1)
    display(['pairConditions : ' num2str(i)])
    temp = allUnitsPSTH_Arr_Base{pairConditions(i,1),ses} - allUnitsPSTH_Arr_Base{pairConditions(i,2),ses}; 
    dataDistanceArr{i,1} = cellfun(@(x) norm(x), mat2cell( temp, size(temp,1), ones(1, size(temp,2)) ));
    clear temp
    temp = allUnitsPSTH_Arr2_Base{pairConditions(i,1),ses} - allUnitsPSTH_Arr2_Base{pairConditions(i,2),ses}; 
    dataDistanceArr2{i,1} = cellfun(@(x) norm(x), mat2cell( temp, size(temp,1), ones(1, size(temp,2)) ));
    clear temp
    ns = 0;
    while ns < numNoise
        ns = ns+1;
        temp = allNoiseUnitPSTH_Arr{pairConditions(i,1),ses,ns} - allNoiseUnitPSTH_Arr{pairConditions(i,2),ses,ns}; 
        noiseDistanceArr{i,ns} = cellfun(@(x) norm(x), mat2cell( temp, size(temp,1), ones(1, size(temp,2)) ));
        clear temp
        temp = allNoiseUnitPSTH_Arr2{pairConditions(i,1),ses,ns} - allNoiseUnitPSTH_Arr2{pairConditions(i,2),ses,ns}; 
        noiseDistanceArr2{i,ns} = cellfun(@(x) norm(x), mat2cell( temp, size(temp,1), ones(1, size(temp,2)) ));
        clear temp
    end
    clear ns
end
clear i

% plot data vs noise
% % % % figure('Name', [Monkey, ', Session ' num2str(ses) ', distance in the high dimentional neural space'])
for timepoint = 1:85
    for pairNum = 1:size(pairConditions,1)
% % % %         subplot(6,6,pairNum), box on
        aa = cell2mat(noiseDistanceArr(pairNum,:)')'; % size(aa) = [numTimeStamps x numNoise]
% % % %         histogram(aa(timepoint,:),'FaceColor','b');
% % % %         line([dataDistanceArr{pairNum}(timepoint),dataDistanceArr{pairNum}(timepoint)],[min(ylim),max(ylim)],'Color',[1,0,0])
% % % %         xlim([0,6])
% % % %         title([{['Time = Cue ' num2str(newT(timepoint))]},{['Distance [EV,Var]: [' num2str(groupingInfo.GroupingParaVal{pairConditions(pairNum,1)}) ']-[' num2str(groupingInfo.GroupingParaVal{pairConditions(pairNum,2)}) ']']}])
% % % %         ylabel('samples')
% % % % % % % %         % wilcoxon
% % % % % % % %         p = ranksum(dataDistanceArr{pairNum}(timepoint),aa(timepoint,:));
% % % % % % % %         % ttest
% % % % % % % %         [~,p] = ttest(aa(timepoint,:), dataDistanceArr{pairNum}(timepoint));
        [p(timepoint, pairNum) , sampleSide] = TestSignificance(aa(timepoint,:), dataDistanceArr{pairNum}(timepoint), 0.05, 'two side');
% % % %         text(0.5, 0.7*max(ylim), ['p = ' num2str(p)])
        clear aa 
    end
% % % %     ax = gca;
% % % %     F(timepoint) = getframe(ax);
end

figure('Name', [Monkey, ', Session ' num2str(ses) ', distance in high dimentional neural space'])
hold on, box on
for i=1:36
    plot( newT, (i*2) + (~isnan(p(:,i))) )
    if ~sum((~isnan(p(:,i)))) 
        txtColor = [1 0 0];
    else
        txtColor = [0 0 0];
    end
    text(-0.09, 2*i, {['\color[rgb]{' num2str(txtColor) '}[' num2str(groupingInfo.GroupingParaVal{pairConditions(i,1)}) ']-[' num2str(groupingInfo.GroupingParaVal{pairConditions(i,2)}) ']']})
end

figure
for i = 1:9
    plot(newT,allUnitsPSTH_Arr_Base{i, 1}(1,:))
    hold on
end
end
%% test trajectories' distance in the low dimentional pc space

if PoolUnitsFromAllSessions
    aa = allUnitsPSTH_Arr_Base;
    for i=1:size(aa,1)
        bb{i,1}=cell2mat(aa(i,:)');
    end
    allUnitsPSTH_Arr_Base = bb;
    clear bb i
end

% apply PCA on this array's data
for ses = 1:size(allUnitsPSTH_Arr_Base,2)
        numTimeStamps = size(allUnitsPSTH_Arr_Base{1,ses},2);
        numUnits =  size(allUnitsPSTH_Arr_Base{1,ses},1);
        obsrv_x_datadimensions = cell2mat(allUnitsPSTH_Arr_Base(:,ses)')'; % numObservations = numGroups x numTimeStamps; numDimensions of each data point = numUnits ; 

        [tempCoeff,tempScore,~,~,tempVar] = pca (obsrv_x_datadimensions);
        pcaSesArr{ses,1}.Coeff = tempCoeff;
        pcaSesArr{ses,1}.Score = tempScore;
        pcaSesArr{ses,1}.PCVar = tempVar;
        pcaSesArr{ses,1}.T = newT;
        pcaSesArr{ses,1}.numTimeStamps = numTimeStamps;
        pcaSesArr{ses,1}.numUnits = numUnits;
        clear obsrv_x_datadimensions numTimeStamps tempCoeff tempScore numUnits tempVar
end
clear ses

if 0 % -----------------------------------------------------
% apply PCA on noise for Array 1
[~ , t2BaseIdx] = min(abs(newT));
allNoiseUnitPSTH_Arr_Base = cellfun(@(x) x - mean(x(:,1:t2BaseIdx),2)  , allNoiseUnitPSTH_Arr, 'UniformOutput', 0);
clear t2BaseIdx
ses = 0;
while ses < size(allNoiseUnitPSTH_Arr_Base,2)
    ses = ses + 1;
    ns = 0;
    while ns < numNoise
        ns = ns + 1;
        numTimeStamps = size(allNoiseUnitPSTH_Arr_Base{1,1},2);
        numUnits      = size(allNoiseUnitPSTH_Arr_Base{1,1},1);
        obsrv_x_datadimensions = cell2mat(allNoiseUnitPSTH_Arr_Base(:,ns)')'; % numObservations = numGroups x numTimeStamps; numDimensions of each data point = numUnits ; 

        [tempCoeff,tempScore] = pca (obsrv_x_datadimensions);
        pcaNoiseSesArr{ses, ns}.Coeff = tempCoeff;
        pcaNoiseSesArr{ses, ns}.Score = tempScore;
        pcaNoiseSesArr{ses, ns}.T = newT;
        pcaNoiseSesArr{ses, ns}.numTimeStamps = numTimeStamps;
        pcaNoiseSesArr{ses, ns}.numUnits = numUnits;
        clear obsrv_x_datadimensions numTimeStamps tempCoeff tempScore numUnits
    end
    clear ns
end
clear ses

% measure distance between conditions for data and noise in the PC space
ses = 11;
% pcsToConsider = [1 2 3];
pcsToConsider = find(cumsum(pcaSesArr{ses}.PCVar)<=75);
for i = 1:size(pairConditions,1)
    display(['pairConditions : ' num2str(i)])
 % Array1  
    % split the data matrix in PC space into conditions: size(conditionsInPC) = [numPCs, numTimeStamps, numConditions]
    conditionsInPC = reshape(pcaSesArr{ses}.Score', size(pcaSesArr{ses}.Score,2), pcaSesArr{ses}.numTimeStamps, size(pcaSesArr{ses}.Score,1)/pcaSesArr{ses}.numTimeStamps);
    % for each time point, subtract data vectors in PC space across two conditions: size(temp) = [numSelectedPCs x numTimeStamps])
    temp = squeeze(conditionsInPC(pcsToConsider, : , pairConditions(i,1))) - squeeze(conditionsInPC(pcsToConsider, : , pairConditions(i,2)));
    % for each time point, get the euclidean distance between two conditions, i.e., length of the subtraction vector
    dataDistanceArr{i,1} = cellfun(@(x) norm(x), mat2cell( temp, size(temp,1), ones(1, size(temp,2)) ));
    clear temp conditionsInPC
 % Array2  
    % split the data matrix in PC space into conditions: size(conditionsInPC) = [numPCs, numTimeStamps, numConditions]
    conditionsInPC = reshape(pcaSesArr2{ses}.Score', size(pcaSesArr2{ses}.Score,2), pcaSesArr2{ses}.numTimeStamps, size(pcaSesArr2{ses}.Score,1)/pcaSesArr2{ses}.numTimeStamps);
    % for each time point, subtract data vectors in PC space across two conditions: size(temp) = [numSelectedPCs x numTimeStamps])
    temp = squeeze(conditionsInPC(pcsToConsider, : , pairConditions(i,1))) - squeeze(conditionsInPC(pcsToConsider, : , pairConditions(i,2)));
    % for each time point, get the euclidean distance between two conditions, i.e., length of the subtraction vector
    dataDistanceArr2{i,1} = cellfun(@(x) norm(x), mat2cell( temp, size(temp,1), ones(1, size(temp,2)) ));
    clear temp conditionsInPC
 % Noise
    ns = 0;
    while ns < numNoise
       ns = ns+1;
       if ~isempty(pcaNoiseSesArr{ses,ns}.Score)
     % Array1  
        % split the data matrix in PC space into conditions: size(conditionsInPC) = [numPCs, numTimeStamps, numConditions]
        conditionsInPC = reshape(pcaNoiseSesArr{ses,ns}.Score', size(pcaNoiseSesArr{ses,ns}.Score,2), pcaNoiseSesArr{ses,ns}.numTimeStamps, size(pcaNoiseSesArr{ses,ns}.Score,1)/pcaNoiseSesArr{ses,ns}.numTimeStamps);
        % for each time point, subtract data vectors in PC space across two conditions: size(temp) = [numSelectedPCs x numTimeStamps])
        temp = squeeze(conditionsInPC(pcsToConsider, : , pairConditions(i,1))) - squeeze(conditionsInPC(pcsToConsider, : , pairConditions(i,2)));
        % for each time point, get the euclidean distance between two conditions, i.e., length of the subtraction vector
        noiseDistanceArr{i,ns} = cellfun(@(x) norm(x), mat2cell( temp, size(temp,1), ones(1, size(temp,2)) ));
        clear temp conditionsInPC
       end
       if ~isempty(pcaNoiseSesArr2{ses,ns}.Score)
     % Array2  
        % split the data matrix in PC space into conditions: size(conditionsInPC) = [numPCs, numTimeStamps, numConditions]
        conditionsInPC = reshape(pcaNoiseSesArr2{ses,ns}.Score', size(pcaNoiseSesArr2{ses,ns}.Score,2), pcaNoiseSesArr2{ses,ns}.numTimeStamps, size(pcaNoiseSesArr2{ses,ns}.Score,1)/pcaNoiseSesArr2{ses,ns}.numTimeStamps);
        % for each time point, subtract data vectors in PC space across two conditions: size(temp) = [numSelectedPCs x numTimeStamps])
        temp = squeeze(conditionsInPC(pcsToConsider, : , pairConditions(i,1))) - squeeze(conditionsInPC(pcsToConsider, : , pairConditions(i,2)));
        % for each time point, get the euclidean distance between two conditions, i.e., length of the subtraction vector
        noiseDistanceArr2{i,ns} = cellfun(@(x) norm(x), mat2cell( temp, size(temp,1), ones(1, size(temp,2)) ));
        clear temp conditionsInPC
       end
    end
    clear ns
end
clear i

% plot data vs noise
% % % % figure('Name', [Monkey, ', Session ' num2str(ses) ', distance in the low dimentional PCs [' num2str(pcsToConsider) '] space'])
for timepoint = 1:85
    for pairNum = 1:size(pairConditions,1)
% % % %         subplot(6,6,pairNum), box on
        aa = cell2mat(noiseDistanceArr(pairNum,:)')'; % size(aa) = [numTimeStamps x numNoise]
% % % %         histogram(aa(timepoint,:),'FaceColor','b');
% % % %         line([dataDistanceArr{pairNum}(timepoint),dataDistanceArr{pairNum}(timepoint)],[min(ylim),max(ylim)],'Color',[1,0,0])
% % % %         xlim([0,6])
% % % %         title([{['Time = Cue ' num2str(newT(timepoint))]},{['Distance [EV,Var]: [' num2str(groupingInfo.GroupingParaVal{pairConditions(pairNum,1)}) ']-[' num2str(groupingInfo.GroupingParaVal{pairConditions(pairNum,2)}) ']']}])
% % % %         ylabel('samples')
% % % % % % % %         % wilcoxon
% % % % % % % %         p = ranksum(dataDistanceArr{pairNum}(timepoint), aa(timepoint,:), 'tail', 'right');
% % % % % % % %         % ttest
% % % % % % % %         [~,p] = ttest(aa(timepoint,:), dataDistanceArr{pairNum}(timepoint));
        [p(timepoint, pairNum) , sampleSide] = TestSignificance(aa(timepoint,:), dataDistanceArr{pairNum}(timepoint), 0.05, 'two side');
% % % %         text(1, 0.7*max(ylim), ['p = ' num2str(p)])
        clear aa 
    end
% % % %     ax = gca;
% % % %     F(timepoint) = getframe(ax);
end
clear timepoint pairNum

figure('Name', [Monkey, ', Session ' num2str(ses) ', distance in the low dimentional PCs [' num2str(pcsToConsider') '] space'])
hold on, box on
for i=1:36
    plot( newT, (i*2) + (~isnan(p(:,i))) )
    if ~sum((~isnan(p(:,i)))) 
        txtColor = [1 0 0];
    else
        txtColor = [0 0 0];
    end
    text(-0.09, 2*i, {['\color[rgb]{' num2str(txtColor) '}[' num2str(groupingInfo.GroupingParaVal{pairConditions(i,1)}) ']-[' num2str(groupingInfo.GroupingParaVal{pairConditions(i,2)}) ']']})
end
clear p

end % --------------------------------------------------------------------
%% plot data in 2D PC space
ses = 1;
figProperties =  struct( 'Name','test', 'Title',[Monkey ', session ' num2str(ses) ', ' groupingCriteria]);
componentsToPlot = [1 2 3];
switch groupingCriteria
    case 'expected_reward'
        groupingInfo.GroupingPara = 'EV';
        colorCode = [{[0 1 1]},...
                     {[0 0 1]},...
                     {[0 0 0]}];
    case 'Variance'
        groupingInfo.GroupingPara = 'Var';
        colorCode = [{[0 1 1]},...
                     {[0 0 1]},...
                     {[0 0 0]}];
    case 'expecter_reward & RewardVariance'
        groupingInfo.GroupingPara = '[EV,Var]';
        colorCode = [{[0.6 1 0.6]},...
                     {[0 1 0]},...
                     {[0.1 0.5 0.1]},...
                     ...
                     {[1 0.6 0.6]},...
                     {[0.99 0.2 0.2]},...
                     {[0.6 0 0]},...
                     ...
                     {[0 1 1]},...
                     {[0 0 1]},...
                     {[0 0 0]},...
                     ];
end
markedTimePoints = struct('Time',{0 0.4 0.8}, 'Marker',[{'o'},{'*'},{'*'}], 'Color',[{'r'},{'g'},{'m'}]);
timeWindowPlot = [0 0.8];

fg = PlotTrajectoriesInStateSpace(pcaSesArr{ses}, '3D', figProperties, componentsToPlot, groupingInfo, markedTimePoints, timeWindowPlot, colorCode);
% -----------------------------------------------

%% plot PSTH for PCs

if PoolUnitsFromAllSessions
    ses = 1;
else 
    ses = 4;
end

componentsToPlot = [1 2 3];
figure, 
legendTxt = [];

for pcNum = 1:length(componentsToPlot)
    subplot(1,length(componentsToPlot),pcNum), hold on, box on
    title(['ExplainedVar = ' num2str(pcaSesArr{ses,1}.PCVar(componentsToPlot(pcNum))) ], 'FontSize', 8)
    xlabel('Time', 'FontWeight', 'Bold')
    if PoolUnitsFromAllSessions
        ylabel([{['PC' num2str(componentsToPlot(pcNum)) '  PSTH']} , {['All Sessions']} ], 'FontWeight', 'Bold')
    else
        ylabel([{['PC' num2str(componentsToPlot(pcNum)) '  PSTH']} , {['Session ' num2str(ses)]} ], 'FontWeight', 'Bold')
    end

    for cnd = 1:size(allUnitsPSTH_Arr_Base,1)
        plot(newT, sum( repmat(pcaSesArr{ses,1}.Coeff(:,componentsToPlot(pcNum)), 1, size(allUnitsPSTH_Arr_Base{cnd,ses},2)) .* allUnitsPSTH_Arr_Base{cnd,ses} ), ...
            'Color', colorCode{1,cnd}, 'LineWidth', 2)
        legendTxt = [legendTxt; {['\color[rgb]{' num2str(colorCode{cnd}) '}' groupingInfo.GroupingPara ' = ' num2str(groupingInfo.GroupingParaVal{cnd})]}];  
    end
end
legend(legendTxt)









