
% 06.05.2018: this is the same scripts as R13_SaveSigUnits_ScatterPlotGLMCoeff but with the new channel-electrode map 

% This scripts applies GLM on spike counts of active units and plots the distribution of the 
% coefficients. also store list of the significant units on the hard drive

% Notice that to run this script, you need the spike count data to be
% already stored on the hard drive. You should use PlotTuningCurve.m to get
% the spike count files and then run CollectAndStoreSpikeCounts

clear, clc
dbstop if error 

spike_dir   = 'X:\data\Add PSTH to Aligned Spk\'; % 'D:\Shared\data\Add PSTH to Aligned Spk'; 
save_dire   = 'X:\Report 16\';  

Monkey = 'Mojo';         % select among 'Mojo', 'MacDuff' 

alignedEvent = '(cue)';     % select among 'cue', 'reward' and    
groupingCriteria = 'CuePos & CurrEV & CurrVar';  % see all possible conditions under GroupTrials()
timeInterval = [0.4,0.8];
minRepetitionPerCondition = 5;
ZScoreEachTuningCurve = 1;

% we need these parameters to read out the right file for list of the active units. these values are included in the file names.
StabilityCheckTimeInterval = [-0.5 1];
AveSpkCountThreshold = 5;

color_code = hsv(9);   % ***** use specific color for each task condition; 36 colors
%% 
% here we use spike count data which is already stored on the hard drive; use CollectAndStoreSpikeCounts.m to get the files in the right format
SpkCountsPath = ['X:\Report 13\Spike Counts\' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']_' groupingCriteria '_minRepetition' num2str(minRepetitionPerCondition) '\'];
% ['D:\Shared\Report 11\Spike Counts\' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']_' groupingCriteria '_minRepetition' num2str(minRepetitionPerCondition) '\'];

% load list of active units for this monkey; this will load variables ActiveSesUnitsArr1 and ActiveSesUnitsArr2: [SesName, ActiveUnitNames, MeanSpikeCount in StabilityCheckTimeInterval]
load(['X:\data\' Monkey '_UnitStability_' alignedEvent '[' regexprep(num2str(StabilityCheckTimeInterval(1)),'0.','0,') ' ' regexprep(num2str(StabilityCheckTimeInterval(2)),'0.','0,') ']_AveSpkCntThreshold' num2str(AveSpkCountThreshold)],...
      'ActiveSesUnitsArr1', 'ActiveSesUnitsArr2');
% get the list of sessions in which there was active units
sesListSpkCnt = unique([ActiveSesUnitsArr1(:,1);ActiveSesUnitsArr2(:,1)]);
%% Apply GLMon individual cells
    
% get the new channel-electrode assignments on the array
eval(['[newArray1, newArray2] = new_electrodepinout' Monkey '(''chan'');']) 

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
    
    % collect mean spike count of all the units, separately for each group 
    if ~exist('allUnitTrlsSpkCountArr1','var')
        allUnitTrlsSpkCountArr1 = []; % each row will be the spike count for one unit
        allUnitTrlsParaValArr1 = [];
        allUnitTrlsUnitIdxArr1 = [];  % indicates each trial sample belongs to which unit and session. indices indicate the unit index in the population not in the session  
        allUnitIDsArray1 = []; % includes session file name and unit lable or each trial
    end
    if ~exist('allUnitTrlsSpkCountArr2','var')
        allUnitTrlsSpkCountArr2 = []; % each row will be the spike count for one unit
        allUnitTrlsParaValArr2 = [];
        allUnitTrlsUnitIdxArr2 = []; % indicates each trial sample belongs to which unit and session. indices indicate the unit index in the population not in the session 
        allUnitIDsArray2 = []; % includes session file name and unit lable or each trial
    end
    
    % collect spike counts from all the units on array 1
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
        
        switch [num2str(thisActiveChanArrayIdx) ' ' num2str(~isempty(spkData.(thisSesActiveUnits{un,1}{1,1})))]
          case '1 1'
            allUnitIDsArray1 = [allUnitIDsArray1; [{sesListSpkCnt{ses}},{thisSesActiveUnits{un,1}{1,1}}]];
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
            
          case '2 1'
          % collect spike counts from all the units on array 2
          % get the list of active units in this session on array 2
            allUnitIDsArray2 = [allUnitIDsArray2; [{sesListSpkCnt{ses}},{thisSesActiveUnits{un,1}{1,1}}]];
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
allUnitTrlsUnitIdxArr1(nanIdx1,:) = [];   
% if units were removed from the spike count, remove it from the unit ID list
allUnitIDsArray1 = allUnitIDsArray1(unique(allUnitTrlsUnitIdxArr1(:,1)),:);

nanIdx2 = isnan(allUnitTrlsSpkCountArr2);
allUnitTrlsSpkCountArr2(nanIdx2) = [];
allUnitTrlsParaValArr2(nanIdx2,:) = [];    
allUnitTrlsUnitIdxArr2(nanIdx2,:) = [];        
% if units were removed from the spike count, remove it from the unit ID list
allUnitIDsArray2 = allUnitIDsArray2(unique(allUnitTrlsUnitIdxArr2(:,1)),:);

clear nanIdx*

% apply GLM on array 1 units ---------------------------------------------------------
uniqUnitIdx1 = unique(allUnitTrlsUnitIdxArr1(:,1));
bGLMUnitArr1 = cell(length(uniqUnitIdx1),1);
pGLMUnitArr1 = cell(length(uniqUnitIdx1),1);
ev = [3 6 9]; zev = zscore(ev);
vr = [0 1 4]; zvr = zscore(vr);
cu = [-1 1];  zcu = zscore(cu);
ev_vr = [0 3 6 9 12 24 36]; z_ev_vr = zscore(ev_vr);

for u = 1:length(uniqUnitIdx1)
    idx = find(uniqUnitIdx1(u) == allUnitTrlsUnitIdxArr1(:,1));
    % zscore expected value for GLM input
    p1_a1 = allUnitTrlsParaValArr1(idx,2);
    for e = 1:length(ev)
        p1_a1(p1_a1 == ev(e)) = zev(e);
    end
    clear e

    % zscore variance for GLM input
    p2_a1 = allUnitTrlsParaValArr1(idx,3);
    for v = 1:length(vr)
        p2_a1(p2_a1 == vr(v)) = zvr(v);
    end
    clear v
    
    % make a new parameter for EV*Var interaction
    p1p2_a1 = allUnitTrlsParaValArr1(idx,2) .* allUnitTrlsParaValArr1(idx,3);
    for evvr = 1:length(ev_vr)
        p1p2_a1(p1p2_a1 == ev_vr(evvr)) = z_ev_vr(evvr);
    end
    clear evvr

    % zscore cue location for GLM input
    p4_a1 = allUnitTrlsParaValArr1(idx,1);
    for c = 1:length(cu)
        p4_a1(p4_a1 == cu(c)) = zcu(c);
    end
    clear c

    %
    [b,~,s] = glmfit([p1_a1 p2_a1 p1p2_a1 p4_a1] , allUnitTrlsSpkCountArr1(idx) , 'normal');
    bGLMUnitArr1{u} = b';
    pGLMUnitArr1{u} = s.p';
    clear s b p1_a1 p2_a1 p1p2_a1 idx p4_a1
end
clear uniqUnitIdx1 u

% apply GLM on array 2 units ---------------------------------------------------------
uniqUnitIdx2 = unique(allUnitTrlsUnitIdxArr2(:,1));
bGLMUnitArr2 = cell(length(uniqUnitIdx2),1);
pGLMUnitArr2 = cell(length(uniqUnitIdx2),1);

for u = 1:length(uniqUnitIdx2)
    idx = find(uniqUnitIdx2(u) == allUnitTrlsUnitIdxArr2(:,1));
    % zscore expected value for GLM input
    p1_a2 = allUnitTrlsParaValArr2(idx,2);
    for e = 1:length(ev)
        p1_a2(p1_a2 == ev(e)) = zev(e);
    end
    clear e

    % zscore variance for GLM input
    p2_a2 = allUnitTrlsParaValArr2(idx,3);
    for v = 1:length(vr)
        p2_a2(p2_a2 == vr(v)) = zvr(v);
    end
    clear v
    
    % make a new parameter for EV*Var interaction
    p1p2_a2 = allUnitTrlsParaValArr2(idx,2) .* allUnitTrlsParaValArr2(idx,3);
    for evvr = 1:length(ev_vr)
        p1p2_a2(p1p2_a2 == ev_vr(evvr)) = z_ev_vr(evvr);
    end
    clear evvr

    % zscore cue location for GLM input
    p4_a2 = allUnitTrlsParaValArr2(idx,1);
    for c = 1:length(cu)
        p4_a2(p4_a2 == cu(c)) = zcu(c);
    end
    clear c

    [b,~,s] = glmfit([p1_a2 p2_a2 p1p2_a2 p4_a2] , allUnitTrlsSpkCountArr2(idx) , 'normal');
    bGLMUnitArr2{u} = b';
    pGLMUnitArr2{u} = s.p';
    clear s b p1_a2 p2_a2 p1p2_a2 idx p4_a2
end
clear uniqUnitIdx2 u

%% Investigate significant units
bGLMUnitArr1 = cell2mat(bGLMUnitArr1);
bGLMUnitArr2 = cell2mat(bGLMUnitArr2);
pGLMUnitArr1 = cell2mat(pGLMUnitArr1);
pGLMUnitArr2 = cell2mat(pGLMUnitArr2);

maxXlim = max([bGLMUnitArr1(:);bGLMUnitArr2(:)]);
minXlim = min([bGLMUnitArr1(:);bGLMUnitArr2(:)]);

% ---- array 1
evSigUnitIdx  = find(pGLMUnitArr1(:,2)<=0.05);
varSigUnitIdx = find(pGLMUnitArr1(:,3)<=0.05);
intSigUnitIdx = find(pGLMUnitArr1(:,4)<=0.05);
cueSigUnitIdx = find(pGLMUnitArr1(:,5)<=0.05);
allSigUnit = unique([evSigUnitIdx; varSigUnitIdx; intSigUnitIdx; cueSigUnitIdx]);
% store significant unit IDs on the hard drive
save([save_dire '\Results\Active Units\' Monkey '_SigUnitIDs_dlPFC_' alignedEvent(2:end-1) ], ...
      'allUnitIDsArray1','evSigUnitIdx','varSigUnitIdx','intSigUnitIdx','cueSigUnitIdx','pGLMUnitArr1','bGLMUnitArr1','timeInterval','alignedEvent',...
      'allUnitTrlsSpkCountArr1', 'allUnitTrlsUnitIdxArr1')

numRows = 2; numColumns = 3;

axisName = [{'EV'}  {'2'};...
            {'Var'} {'3'};...
            {'Int'} {'4'};...
            {'Cue'} {'5'}];

axisPairs = nchoosek(1:size(axisName,1),2);       


fg1 = figure('Name', [Monkey ', GLM significant units'], 'Position', [20 60 0 0] + (0.8*get(0,'ScreenSize')));
for ax = 1:size(axisPairs,1)
    subplot(numRows,numColumns,ax), hold on, box on   
    xlabel([axisName{axisPairs(ax,1),1} ' coefficient'], 'FontWeight', 'bold', 'FontSize', 14)
    ylabel([axisName{axisPairs(ax,2),1} ' coefficient'], 'FontWeight', 'bold', 'FontSize', 14)
    plot([-2 2],[-2 2],'Color', [0.7 0.7 0.7]), axis square, title(['dlPFC, ' alignedEvent ' [' num2str(timeInterval) ']'])
    set(gca, 'XLim' , [minXlim,maxXlim], 'YLim' , [minXlim,maxXlim])

    eval(['scatter( bGLMUnitArr1(evSigUnitIdx,'  axisName{axisPairs(ax,1),2} ') , bGLMUnitArr1(evSigUnitIdx,' axisName{axisPairs(ax,2),2} ')  , ''o'', ''MarkerEdgeColor'', [1 0 0] , ''LineWidth'',1.5); '])
    eval(['scatter( bGLMUnitArr1(varSigUnitIdx,' axisName{axisPairs(ax,1),2} ') , bGLMUnitArr1(varSigUnitIdx,' axisName{axisPairs(ax,2),2} ') , 20, ''*'', ''MarkerEdgeColor'', [0 0.7 0.7] );'])
    eval(['scatter( bGLMUnitArr1(intSigUnitIdx,' axisName{axisPairs(ax,1),2} ') , bGLMUnitArr1(intSigUnitIdx,' axisName{axisPairs(ax,2),2} ') , 140, ''s'', ''MarkerEdgeColor'', [0 0 0] );'])
    eval(['scatter( bGLMUnitArr1(cueSigUnitIdx,' axisName{axisPairs(ax,1),2} ') , bGLMUnitArr1(cueSigUnitIdx,' axisName{axisPairs(ax,2),2} ') , 100, ''d'', ''MarkerEdgeColor'', ''m'' );'])
    legend('unity line','EV significant', 'Var significant','EV*Var significant','CueLocation signi.')

end

saveas(fg1, [save_dire  'Results\Active Units\' Monkey '_SignificantUnits_GLMcoeffScatter_dlPFC'], 'fig')
saveas(fg1, [save_dire  'Results\Active Units\' Monkey '_SignificantUnits_GLMcoeffScatter_dlPFC'], 'bmp')
close(fg1)
clear evSigUnitIdx varSigUnitIdx intSigUnitIdx allSigUnit cueSigUnitIdx

% ---- array 2
evSigUnitIdx  = find(pGLMUnitArr2(:,2)<=0.05);
varSigUnitIdx = find(pGLMUnitArr2(:,3)<=0.05);
intSigUnitIdx = find(pGLMUnitArr2(:,4)<=0.05);
cueSigUnitIdx = find(pGLMUnitArr2(:,5)<=0.05);
allSigUnit = unique([evSigUnitIdx; varSigUnitIdx; intSigUnitIdx; cueSigUnitIdx]);
% store significant unit IDs on the hard drive
save([save_dire '\Results\Active Units\' Monkey '_SigUnitIDs_7a_' alignedEvent(2:end-1) ], ...
      'allUnitIDsArray2','evSigUnitIdx','varSigUnitIdx','intSigUnitIdx','cueSigUnitIdx','pGLMUnitArr2','bGLMUnitArr2','timeInterval','alignedEvent',...
      'allUnitTrlsSpkCountArr2', 'allUnitTrlsUnitIdxArr2')


fg2 = figure('Name', [Monkey ', GLM significant units'], 'Position', [20 60 0 0] + (0.8*get(0,'ScreenSize')));
for ax = 1:size(axisPairs,1)
    subplot(numRows,numColumns,ax), hold on, box on   
    xlabel([axisName{axisPairs(ax,1),1} ' coefficient'], 'FontWeight', 'bold', 'FontSize', 14)
    ylabel([axisName{axisPairs(ax,2),1} ' coefficient'], 'FontWeight', 'bold', 'FontSize', 14)
    plot([-2 2],[-2 2],'Color', [0.7 0.7 0.7]), axis square, title(['7a, ' alignedEvent ' [' num2str(timeInterval) ']'])
    set(gca, 'XLim' , [minXlim,maxXlim], 'YLim' , [minXlim,maxXlim])

    eval(['scatter( bGLMUnitArr2(evSigUnitIdx,'  axisName{axisPairs(ax,1),2} ') , bGLMUnitArr2(evSigUnitIdx,' axisName{axisPairs(ax,2),2} ')  , ''o'', ''MarkerEdgeColor'', [1 0 0] , ''LineWidth'',1.5); '])
    eval(['scatter( bGLMUnitArr2(varSigUnitIdx,' axisName{axisPairs(ax,1),2} ') , bGLMUnitArr2(varSigUnitIdx,' axisName{axisPairs(ax,2),2} ') , 20, ''*'', ''MarkerEdgeColor'', [0 0.7 0.7] );'])
    eval(['scatter( bGLMUnitArr2(intSigUnitIdx,' axisName{axisPairs(ax,1),2} ') , bGLMUnitArr2(intSigUnitIdx,' axisName{axisPairs(ax,2),2} ') , 140, ''s'', ''MarkerEdgeColor'', [0 0 0] );'])
    eval(['scatter( bGLMUnitArr2(cueSigUnitIdx,' axisName{axisPairs(ax,1),2} ') , bGLMUnitArr2(cueSigUnitIdx,' axisName{axisPairs(ax,2),2} ') , 100, ''d'', ''MarkerEdgeColor'', ''m'' );'])
    legend('unity line','EV significant', 'Var significant','EV*Var significant','CueLocation signi.')

end

saveas(fg2, [save_dire  'Results\Active Units\' Monkey '_SignificantUnits_GLMcoeffScatter_7a'], 'fig')
saveas(fg2, [save_dire  'Results\Active Units\' Monkey '_SignificantUnits_GLMcoeffScatter_7a'], 'bmp')
close(fg2)
clear evSigUnitIdx varSigUnitIdx intSigUnitIdx allSigUnit cueSigUnitIdx

