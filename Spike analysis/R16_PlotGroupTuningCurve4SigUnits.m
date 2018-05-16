
% 09.05.2018: this is the same scripts as R14_PlotGroupTuningCurve4SigUnits but with the new channel-electrode map 

% this script plots and stores tuning curve of GLM-significant active units for both areas and all GLM factors to 9 [EV Var] 
% task conditions.
% To run this script we need to know which of the active units had significant GLM cofficients. To get this you need to 
% run: 1- CollectAndStoreSpikeCounts and 2- R13_SaveSigUnits_ScatterPlotGLMCoeff.m


clear, clc
dbstop if error 

Monkey = 'Mojo';         % select among 'Mojo', 'MacDuff' 

alignedEvent = '(cue)';     % select among 'cue', 'reward' and    
groupingCriteria = 'CuePos & CurrEV & CurrVar';  % see all possible conditions under GroupTrials()
timeInterval = [0.2,0.8];
minRepetitionPerCondition = 5;
ZScoreEachTuningCurve = 1;

colorCode = [[0 1 1]; ...
             [0 0 1]; ...
             [0 0 0]];

SpkCountsPath = ['X:\Report 13\Spike Counts\' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']_' groupingCriteria '_minRepetition' num2str(minRepetitionPerCondition) '\'];
SignificantUnitsPath = 'X:\Report 16';
ParaNames = {'ev', 'var', 'int'};
Area = {'dlPFC','7a'};

numRows = 2;
numColumns = 3;

save_fig_dir = 'X:\Report 16\Grouped Tuning Curve\';
%%
% get the new channel-electrode assignments on the array
eval(['[newArray1, newArray2] = new_electrodepinout' Monkey '(''chan'');']) 

for areA = 1:length(Area)
    % load Significant Unit IDs for Area; this will load following variables:
    %       'allUnitIDsArray1','evSigUnitIdx','varSigUnitIdx','intSigUnitIdx','cueSigUnitIdx','pGLMUnitArr1','bGLMUnitArr1','timeInterval','alignedEvent','allUnitTrlsSpkCountArr1', 'allUnitTrlsUnitIdxArr1'
    %  allUnitTrlsUnitIdxArr1: [unit?IdxInPop sesIdx]
    load([SignificantUnitsPath '\Results\Active Units\' Monkey '_SigUnitIDs_' Area{areA} '_' alignedEvent(2:end-1)])
    
    % see what parameters each unit is significant for
    eval(['allUnitsSignificance = nan(size(allUnitIDsArray' num2str(areA) ',1), 3);'])
    
    eval(['allUnitsSignificance(evSigUnitIdx,1)   = sign( bGLMUnitArr' num2str(areA) '(evSigUnitIdx, 2) );']);  % -> [EV  Var  Int  CueLocation]
    eval(['allUnitsSignificance(varSigUnitIdx,2)  = sign( bGLMUnitArr' num2str(areA) '(varSigUnitIdx,3) );']);
    eval(['allUnitsSignificance(intSigUnitIdx,3)  = sign( bGLMUnitArr' num2str(areA) '(intSigUnitIdx,4) );']);
    
    eval(['bGLM = bGLMUnitArr' num2str(areA) ';'])  % -> [Constant  EV  Var  Int  CueLocation]
    unsOfInterest{1} = find( bGLM(:,2)>=0  &  bGLM(:,3)>=0  &  bGLM(:,4)<=0  &  ~isnan(allUnitsSignificance(:,3)) ); 
    unsOfInterestInfo{1} = ['signi.Interaction & bEV>=0 & bVar>=0 & bInt<=0'];

    unsOfInterest{2} = find( bGLM(:,2)>=0  &  bGLM(:,3)<=0  &  bGLM(:,4)<=0  &  ~isnan(allUnitsSignificance(:,3)) ); 
    unsOfInterestInfo{2} = ['signi.Interaction & bEV>=0 & bVar<=0 & bInt<=0'];
    
    unsOfInterest{3} = find( bGLM(:,2)<=0  &  bGLM(:,3)>=0  &  bGLM(:,4)>=0  &  ~isnan(allUnitsSignificance(:,3)) ); 
    unsOfInterestInfo{3} = ['signi.Interaction & bEV<=0 & bVar>=0 & bInt>=0'];
    
    unsOfInterest{4} = find( bGLM(:,2)<=0  &  bGLM(:,3)<=0  &  bGLM(:,4)>=0  &  ~isnan(allUnitsSignificance(:,3)) ); 
    unsOfInterestInfo{4} = ['signi.Interaction & bEV<=0 & bVar<=0 & bInt>=0'];
    
    unsOfInterestOutput = cell(1,length(unsOfInterest));

    for unitType = 1:length(unsOfInterest)
        for un = 1:length(unsOfInterest{unitType})
            display(['unitType: ' num2str(unitType) '   unit: ' num2str(un)])
            idx = unsOfInterest{unitType}(un);
% % % %             eval(['load([SpkCountsPath allUnitIDsArray' num2str(areA) '{idx,1}] , [''spkDataArr'' num2str(areA) ])'])
            eval(['load([SpkCountsPath allUnitIDsArray' num2str(areA) '{idx,1}] , ''spkDataArr1'' , ''spkDataArr2'')'])
            % merge spkDataArr1 and spkDataArr2 fields into one structure
            mrgSpkData = cell2struct([struct2cell(spkDataArr1); struct2cell(spkDataArr2)],  [fieldnames(spkDataArr1);fieldnames(spkDataArr2)], 1);

            eval(['output = mrgSpkData.(allUnitIDsArray' num2str(areA) '{idx,2});'])

            % convert "CuePos & CurrEV & CurrVar" grouping to "CurrEV & CurrVar" grouping
            % and make newOutput with 9 conditions 
            for newGrp = 1:round(length(output)/2)
                newOutput(newGrp).TrialIdx = cell2mat({output(2*newGrp-1:2*newGrp).TrialIdx}');
                newOutput(newGrp).Value = output(2*newGrp-1).Value(2:3);
                newOutput(newGrp).GroupingType = output(2*newGrp-1).GroupingType(10:end);
                newOutput(newGrp).UnitActiveTrials = cell2mat({output(2*newGrp-1:2*newGrp).UnitActiveTrials}');
                newOutput(newGrp).RepPerCnd = output(2*newGrp-1).RepPerCnd + output(2*newGrp).RepPerCnd;
                newOutput(newGrp).AlignedSpikeTimes = [output(2*newGrp-1).AlignedSpikeTimes;output(2*newGrp).AlignedSpikeTimes];
                newOutput(newGrp).SpikeCount = cell2mat({output(2*newGrp-1:2*newGrp).SpikeCount}');
            end
            output = newOutput;
            clear newOutput newGrp
            
            unsOfInterestOutput{unitType} = [unsOfInterestOutput{unitType}; output'];
            clear output idx
        end 
        clear un mrgSpkData spkDataArr*
    end
    clear unitType
    
    
    for unitType = 1:length(unsOfInterest)
        if isempty(unsOfInterestOutput{unitType})
            continue
        end
        fg = figure; hold on, box on
        xlabel('EV', 'FontWeight','Bold'); 
        ylabel('Mean spike count', 'FontWeight','Bold')
        title([ {[Monkey ',  ' 'Area ' Area{areA} ', ' num2str(length(unsOfInterest{unitType})) ' Units,  ' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']']},...
                {unsOfInterestInfo{unitType}}], 'FontWeight','Normal', 'FontSize',8)

        % get the parameter values for each group
        [paraVal,ia,ic] = unique( cell2mat({unsOfInterestOutput{unitType}.Value}') , 'rows' );
        aa = cell2mat({unsOfInterestOutput{unitType}.Value}');
        txt = [];
        for vr = [0 1 4]
            % initialize cuvVarX
            eval(['crvVar' num2str(vr) ' = [];']);
            eval(['semVar' num2str(vr) ' = [];']);
            
            % find the group parameters with this variance value
            idxVar = find(paraVal(:,2) == vr);
            for eevv = 1:length(idxVar) % these are Ev 3, 6, 9
                bb = find(ismember(aa, [paraVal(idxVar(eevv),1) vr], 'rows'));
                eval([ 'crvVar' num2str(vr) ' = [crvVar' num2str(vr) ', '  num2str(mean( cell2mat({unsOfInterestOutput{1,unitType}(bb).SpikeCount}') )) '];' ]);
                eval([ 'semVar' num2str(vr) ' = [semVar' num2str(vr) ', '  num2str(std( cell2mat({unsOfInterestOutput{1,unitType}(bb).SpikeCount}') )/length( cell2mat({unsOfInterestOutput{1,unitType}(bb).SpikeCount}') )) '];' ]);    
                clear bb
            end
            clear idxVar eevv
            
            eval([ 'errorbar( 1:3, crvVar' num2str(vr) ', semVar' num2str(vr) ', ''Color'', colorCode(sqrt(vr)+1,:) )'])
            txt = [txt; {['\color[rgb]{' num2str(colorCode(sqrt(vr)+1, :)) '} Var = ' num2str(vr)]}];
        end
        
        set(gca, 'XTick', 1:3, 'XTickLabel', [{'3'},{'6'},{'9'}], 'XLim', [0 4])
        clear sem* crv* paraVal txt vr bGLM* aa ia ic 
        eval(['clear spkDataArr' num2str(areA) ])
        
        saveas(fg, [save_fig_dir Monkey '_UnitsOfType' num2str(unitType) '_' Area{areA}], 'fig')
        saveas(fg, [save_fig_dir Monkey '_UnitsOfType' num2str(unitType) '_' Area{areA}], 'bmp')
        close(fg); 
        clear fg
    end % end of unitType = 1:length(unsOfInterest)
    clear allUnitsSignificance unitType bGLM unsOfInterest* allUnitIDsArray* *SigUnitIdx
    
end % end of areA
    
    
    