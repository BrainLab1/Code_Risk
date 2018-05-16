
% This scripts applies GLM on spike counts of active units and plots the distribution of the 
% coefficients. 

% 06.05.2018: this is the same scripts as R13_ApplyGLMonSpikeCount but with the new channel-electrode map 
% 1- Notice that to run this script, you need the spike count data to be
% already stored on the hard drive. You should use PlotTuningCurve.m to get
% the spike count files and then run CollectAndStoreSpikeCounts
% 2- The only difference between R13_ApplyGLM... and R12_ApplyGLM... is the
% in R13_ we included CueLocation as a factor in the GLM

clear, clc
dbstop if error 

spike_dir   = 'X:\data\Add PSTH to Aligned Spk\'; % 'Spike_PSTH_ValueTask\Spike_PSTH_CondCodeAdded_ValueTask\';

Monkey = 'Mojo';         % select among 'Mojo', 'MacDuff' 

alignedEvent = '(cue)';     % select among 'cue', 'reward' and    
groupingCriteria = 'CuePos & CurrEV & CurrVar';  % see all possible conditions under GroupTrials()
timeInterval = [0.4,0.8];
minRepetitionPerCondition = 5;
ZScoreEachTuningCurve = 1;

color_code = hsv(9);   % ***** use specific color for each task condition; 36 colors

ApplyGLMonPopulation = 0;
ApplyGLMonIndividualUnits = 1 ;

StabilityCheckTimeInterval = [-0.5 1];
AveSpkCountThreshold = 5;

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

%%
if ApplyGLMonIndividualUnits

% here we use spike count data which is already stored on the hard drive
SpkCountsPath = ['X:\Report 13\Spike Counts\' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']_' groupingCriteria '_minRepetition' num2str(minRepetitionPerCondition) '\'];

% load list of active units for this monkey
load(['X:\data\' Monkey '_UnitStability_' alignedEvent '[' regexprep(num2str(StabilityCheckTimeInterval(1)),'0.','0,') ' ' regexprep(num2str(StabilityCheckTimeInterval(2)),'0.','0,') ']_AveSpkCntThreshold' num2str(AveSpkCountThreshold)],...
      'ActiveSesUnitsArr1', 'ActiveSesUnitsArr2');

sesListSpkCnt = unique([ActiveSesUnitsArr1(:,1);ActiveSesUnitsArr2(:,1)]);

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
        allUnitTrlsUnitIdxArr1 = [];
    end
    if ~exist('allUnitTrlsSpkCountArr2','var')
        allUnitTrlsSpkCountArr2 = []; % each row will be the spike count for one unit
        allUnitTrlsParaValArr2 = [];
        allUnitTrlsUnitIdxArr2 = [];
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
allUnitTrlsUnitIdxArr1(nanIdx1) = [];        
    
nanIdx2 = isnan(allUnitTrlsSpkCountArr2);
allUnitTrlsSpkCountArr2(nanIdx2) = [];
allUnitTrlsParaValArr2(nanIdx2,:) = [];    
allUnitTrlsUnitIdxArr2(nanIdx2) = [];        

clear nanIdx*

% apply GLM on array 1 units ---------------------------------------------------------
% initialize variables
uniqUnitIdx1 = unique(allUnitTrlsUnitIdxArr1);
bGLMUnitArr1 = cell(length(uniqUnitIdx1),1);
pGLMUnitArr1 = cell(length(uniqUnitIdx1),1);
ev = [3 6 9]; zev = zscore(ev);
vr = [0 1 4]; zvr = zscore(vr);
cu = [-1 1];  zcu = zscore(cu);
ev_vr = [0 3 6 9 12 24 36]; z_ev_vr = zscore(ev_vr);

for u = 1:length(uniqUnitIdx1)
    idx = find(uniqUnitIdx1(u) == allUnitTrlsUnitIdxArr1);
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
% initialize variables
uniqUnitIdx2 = unique(allUnitTrlsUnitIdxArr2);
bGLMUnitArr2 = cell(length(uniqUnitIdx2),1);
pGLMUnitArr2 = cell(length(uniqUnitIdx2),1);

for u = 1:length(uniqUnitIdx2)
    idx = find(uniqUnitIdx2(u) == allUnitTrlsUnitIdxArr2);
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

    %
    [b,~,s] = glmfit([p1_a2 p2_a2 p1p2_a2 p4_a2] , allUnitTrlsSpkCountArr2(idx) , 'normal');
    bGLMUnitArr2{u} = b';
    pGLMUnitArr2{u} = s.p';
    clear s b p1_a2 p2_a2 p1p2_a2 idx p4_a2
end
clear uniqUnitIdx2 u

% plot GLM resutls
bGLMUnitArr1 = cell2mat(bGLMUnitArr1);
bGLMUnitArr2 = cell2mat(bGLMUnitArr2);
pGLMUnitArr1 = cell2mat(pGLMUnitArr1);
pGLMUnitArr2 = cell2mat(pGLMUnitArr2);

% % % % binLeftEdge = (-0.5:0.01:0.5-0.005)-0.005;
maxXlim = max([bGLMUnitArr1(:);bGLMUnitArr2(:)]);
minXlim = min([bGLMUnitArr1(:);bGLMUnitArr2(:)]);
binLeftEdge = (minXlim:0.01:maxXlim-0.005)-0.005;
% get the histogram of EV coefficients
nEVArr1 = histcounts( bGLMUnitArr1(:,2), 'BinEdges', binLeftEdge);
nEVArr2 = histcounts( bGLMUnitArr2(:,2), 'BinEdges', binLeftEdge);
% get the histogram of Var coefficients
nVarArr1 = histcounts( bGLMUnitArr1(:,3), 'BinEdges', binLeftEdge);
nVarArr2 = histcounts( bGLMUnitArr2(:,3), 'BinEdges', binLeftEdge);
%
nIntArr1 = histcounts( bGLMUnitArr1(:,4), 'BinEdges', binLeftEdge);
nIntArr2 = histcounts( bGLMUnitArr2(:,4), 'BinEdges', binLeftEdge);

nCueArr1 = histcounts( bGLMUnitArr1(:,5), 'BinEdges', binLeftEdge);
nCueArr2 = histcounts( bGLMUnitArr2(:,5), 'BinEdges', binLeftEdge);

% get the histogram of significant EV coefficients
nSigEVArr1 = histcounts( bGLMUnitArr1(pGLMUnitArr1(:,2)<=0.05,2), 'BinEdges', binLeftEdge);
nSigEVArr2 = histcounts( bGLMUnitArr2(pGLMUnitArr2(:,2)<=0.05,2), 'BinEdges', binLeftEdge);
% get the histogram of significant Var coefficients
nSigVarArr1 = histcounts( bGLMUnitArr1(pGLMUnitArr1(:,3)<=0.05,3), 'BinEdges', binLeftEdge);
nSigVarArr2 = histcounts( bGLMUnitArr2(pGLMUnitArr2(:,3)<=0.05,3), 'BinEdges', binLeftEdge);
%
nSigIntArr1 = histcounts( bGLMUnitArr1(pGLMUnitArr1(:,4)<=0.05,4), 'BinEdges', binLeftEdge);
nSigIntArr2 = histcounts( bGLMUnitArr2(pGLMUnitArr2(:,4)<=0.05,4), 'BinEdges', binLeftEdge);
%
nSigCueArr1 = histcounts( bGLMUnitArr1(pGLMUnitArr1(:,5)<=0.05,5), 'BinEdges', binLeftEdge);
nSigCueArr2 = histcounts( bGLMUnitArr2(pGLMUnitArr2(:,5)<=0.05,5), 'BinEdges', binLeftEdge);


txtLocation = [-0.4, 4]; yAxisLim = [0 7];
fg = figure('Name', [Monkey ', GLM on individual units'], 'Position', [1965         138        1859         787]);
numRows = 2; numClmns = 4;
text(-0.1,0.1,['\angle GLM on dlPFC units, ' alignedEvent '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']'])
subplot(numRows,numClmns,1), title(['EV coefficient'], 'FontSize', 16)
hold on, box on
xlabel('Coeff. value')
ylabel('% Units')
% plot histogram of EV coefficients on array 1
bar(binLeftEdge(1:end-1)+0.005, 100*nEVArr1/sum(nEVArr1), 'FaceColor', [200 200 200]/255, 'EdgeColor', 'none')
bar(binLeftEdge(1:end-1)+0.005, 100*nSigEVArr1/sum(nEVArr1), 'FaceColor', 'b', 'EdgeColor', 'none')
set(gca, 'XLim', [minXlim,maxXlim], 'YLim', yAxisLim)
plot([0 0],ylim,'k')
text(-0.7,0.01,['GLM on dlPFC units, ' alignedEvent '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']'],'Rotation',90, 'FontWeight', 'bold')
text(txtLocation(1), txtLocation(2), [ {['N = ' num2str(size(bGLMUnitArr1,1))]},...
     {['Nsignificant = ' num2str(sum(pGLMUnitArr1(:,2)<=0.05))]},{['%significant = ' num2str(100*sum(pGLMUnitArr1(:,2)<=0.05)/size(bGLMUnitArr1,1))]} ], ...
     'FontSize', 14)
% plot histogram of Var coefficients on array 1
subplot(numRows,numClmns,2), hold on, box on, title('Var coefficient', 'FontSize', 16)
xlabel('Coeff. value')
ylabel('% Units')
bar(binLeftEdge(1:end-1)+0.005, 100*nVarArr1/sum(nVarArr1), 'FaceColor', [200 200 200]/255, 'EdgeColor', 'none')
bar(binLeftEdge(1:end-1)+0.005, 100*nSigVarArr1/sum(nVarArr1), 'FaceColor', 'm', 'EdgeColor', 'none')
set(gca, 'XLim', [minXlim,maxXlim], 'YLim', yAxisLim)
plot([0 0],ylim,'k')
text(txtLocation(1), txtLocation(2), [ {['N = ' num2str(size(bGLMUnitArr1,1))]},...
     {['Nsignificant = ' num2str(sum(pGLMUnitArr1(:,3)<=0.05))]},{['%significant = ' num2str(100*sum(pGLMUnitArr1(:,3)<=0.05)/size(bGLMUnitArr1,1))]} ], ...
     'FontSize', 14)
% plot histogram of EV*Var coefficients on array 1
subplot(numRows,numClmns,3), hold on, box on, title('EV*Var coefficient', 'FontSize', 16)
xlabel('Coeff. value')
ylabel('% Units')
bar(binLeftEdge(1:end-1)+0.005, 100*nIntArr1/sum(nIntArr1), 'FaceColor', [200 200 200]/255, 'EdgeColor', 'none')
bar(binLeftEdge(1:end-1)+0.005, 100*nSigIntArr1/sum(nIntArr1), 'FaceColor', 'r', 'EdgeColor', 'none')
set(gca, 'XLim', [minXlim,maxXlim], 'YLim', yAxisLim)
plot([0 0],ylim,'k')
text(txtLocation(1), txtLocation(2), [ {['N = ' num2str(size(bGLMUnitArr1,1))]},...
     {['Nsignificant = ' num2str(sum(pGLMUnitArr1(:,4)<=0.05))]},{['%significant = ' num2str(100*sum(pGLMUnitArr1(:,4)<=0.05)/size(bGLMUnitArr1,1))]} ], ...
     'FontSize', 14)
% plot histogram of Cue Location coefficients on array 1
subplot(numRows,numClmns,4), hold on, box on, title('Cue Location coefficient', 'FontSize', 16)
xlabel('Coeff. value')
ylabel('% Units')
bar(binLeftEdge(1:end-1)+0.005, 100*nCueArr1/sum(nCueArr1), 'FaceColor', [200 200 200]/255, 'EdgeColor', 'none')
bar(binLeftEdge(1:end-1)+0.005, 100*nSigCueArr1/sum(nCueArr1), 'FaceColor', 'g', 'EdgeColor', 'none')
set(gca, 'XLim', [minXlim,maxXlim], 'YLim', yAxisLim)
plot([0 0],ylim,'k')
text(txtLocation(1), txtLocation(2), [ {['N = ' num2str(size(bGLMUnitArr1,1))]},...
     {['Nsignificant = ' num2str(sum(pGLMUnitArr1(:,5)<=0.05))]},{['%significant = ' num2str(100*sum(pGLMUnitArr1(:,5)<=0.05)/size(bGLMUnitArr1,1))]} ], ...
     'FontSize', 14)



subplot(numRows,numClmns,5), title('EV coefficient', 'FontSize', 16)
hold on, box on
% plot histogram of EV coefficients on array 2
xlabel('Coeff. value')
ylabel('% Units')
bar(binLeftEdge(1:end-1)+0.005, 100*nEVArr2/sum(nEVArr2), 'FaceColor', [200 200 200]/255, 'EdgeColor', 'none')
bar(binLeftEdge(1:end-1)+0.005, 100*nSigEVArr2/sum(nEVArr2), 'FaceColor', 'b', 'EdgeColor', 'none')
set(gca, 'XLim', [minXlim,maxXlim], 'YLim', yAxisLim)
plot([0 0],ylim,'k')
text(-0.7,0.01,['GLM on 7a units, ' alignedEvent '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']'],'Rotation',90, 'FontWeight', 'bold')
text(txtLocation(1), txtLocation(2), [ {['N = ' num2str(size(bGLMUnitArr2,1))]},...
     {['Nsignificant = ' num2str(sum(pGLMUnitArr2(:,2)<=0.05))]},{['%significant = ' num2str(100*sum(pGLMUnitArr2(:,2)<=0.05)/size(bGLMUnitArr2,1))]} ], ...
     'FontSize', 14)
% plot histogram of Var coefficients on array 2
subplot(numRows,numClmns,6), hold on, box on, title('Var coefficient', 'FontSize', 16)
xlabel('Coeff. value')
ylabel('% Units')
bar(binLeftEdge(1:end-1)+0.005, 100*nVarArr2/sum(nVarArr2), 'FaceColor', [200 200 200]/255, 'EdgeColor', 'none')
bar(binLeftEdge(1:end-1)+0.005, 100*nSigVarArr2/sum(nVarArr2), 'FaceColor', 'm', 'EdgeColor', 'none')
set(gca, 'XLim', [minXlim,maxXlim], 'YLim', yAxisLim)
plot([0 0],ylim,'k')
text(txtLocation(1), txtLocation(2), [ {['N = ' num2str(size(bGLMUnitArr2,1))]},...
     {['Nsignificant = ' num2str(sum(pGLMUnitArr2(:,3)<=0.05))]},{['%significant = ' num2str(100*sum(pGLMUnitArr2(:,3)<=0.05)/size(bGLMUnitArr2,1))]} ], ...
     'FontSize', 14)
% plot histogram of EV*Var coefficients on array 2
subplot(numRows,numClmns,7), hold on, box on, title('EV*Var coefficient', 'FontSize', 16)
xlabel('Coeff. value')
ylabel('% Units')
bar(binLeftEdge(1:end-1)+0.005, 100*nIntArr2/sum(nIntArr2), 'FaceColor', [200 200 200]/255, 'EdgeColor', 'none')
bar(binLeftEdge(1:end-1)+0.005, 100*nSigIntArr2/sum(nIntArr2), 'FaceColor', 'r', 'EdgeColor', 'none')
set(gca, 'XLim', [minXlim,maxXlim], 'YLim', yAxisLim)
plot([0 0],ylim,'k')
text(txtLocation(1), txtLocation(2), [ {['N = ' num2str(size(bGLMUnitArr2,1))]},...
     {['Nsignificant = ' num2str(sum(pGLMUnitArr2(:,4)<=0.05))]},{['%significant = ' num2str(100*sum(pGLMUnitArr2(:,4)<=0.05)/size(bGLMUnitArr2,1))]} ], ...
     'FontSize', 14)
% plot histogram of Cue Location coefficients on array 2
subplot(numRows,numClmns,8), hold on, box on, title('Cue Location coefficient', 'FontSize', 16)
xlabel('Coeff. value')
ylabel('% Units')
bar(binLeftEdge(1:end-1)+0.005, 100*nCueArr2/sum(nCueArr2), 'FaceColor', [200 200 200]/255, 'EdgeColor', 'none')
bar(binLeftEdge(1:end-1)+0.005, 100*nSigCueArr2/sum(nCueArr2), 'FaceColor', 'g', 'EdgeColor', 'none')
set(gca, 'XLim', [minXlim,maxXlim], 'YLim', yAxisLim)
plot([0 0],ylim,'k')
text(txtLocation(1), txtLocation(2), [ {['N = ' num2str(size(bGLMUnitArr2,1))]},...
     {['Nsignificant = ' num2str(sum(pGLMUnitArr2(:,5)<=0.05))]},{['%significant = ' num2str(100*sum(pGLMUnitArr2(:,5)<=0.05)/size(bGLMUnitArr2,1))]} ], ...
     'FontSize', 14)


saveas(fg, [SpkCountsPath  'Results\Active Units\' Monkey '_GLMcoeffDist_' groupingCriteria], 'fig')
saveas(fg, [SpkCountsPath  'Results\Active Units\' Monkey '_GLMcoeffDist_' groupingCriteria], 'bmp')
set(gcf, 'Renderer', 'painter')
saveas(fg, [SpkCountsPath  'Results\Active Units\' Monkey '_GLMcoeffDist_' groupingCriteria], 'eps')

end % end of ApplyGLMonIndividualUnits
