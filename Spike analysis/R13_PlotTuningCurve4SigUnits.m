
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
ParaNames = {'ev', 'var', 'int'};
Area = {'dlPFC','7a'};

numRows = 2;
numColumns = 3;

%%
for areA = 1:length(Area)
    % load Significant Unit IDs for Area; this will load following variables:
    %       'allUnitIDsArray1','evSigUnitIdx','varSigUnitIdx','intSigUnitIdx','cueSigUnitIdx','pGLMUnitArr1','bGLMUnitArr1','timeInterval','alignedEvent','allUnitTrlsSpkCountArr1', 'allUnitTrlsUnitIdxArr1'
    %  allUnitTrlsUnitIdxArr1: [unit?IdxInPop sesIdx]
    load([SpkCountsPath '\Results\Active Units\' Monkey '_SigUnitIDs_' Area{areA} '_' alignedEvent(2:end-1)])

    sigPara = 0;
    while sigPara<length(ParaNames) % for each GLM factor
        sigPara = sigPara+1;
        % load session and unit IDs which had significant GLM coeff. for this factor
        eval(['sigSessions = unique(allUnitIDsArray' num2str(areA) '(' ParaNames{sigPara} 'SigUnitIdx,1));'])
        eval(['sigSesAndUnits = allUnitIDsArray' num2str(areA) '(' ParaNames{sigPara} 'SigUnitIdx,:);'])
        
        % initialize figure
        fgCount = 1;
        fg = figure('Name',[Monkey '_' upper(ParaNames{sigPara}) ' Significant Units_' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']'], 'Position', [40 70 0 0] + 0.8*get(0,'ScreenSize'));
        subPltCount = 0;

        for ses = 1:length(sigSessions)
            % load spkDataArr' num2str(areA) ' for this significant sessions; this will load grouped spike counts and task parameter values for active units. 
            load([SpkCountsPath sigSessions{ses}] , ['spkDataArr' num2str(areA) ])
            % read out list of significant units in this sessoin
            thisSesSigUnits =  sigSesAndUnits( find(strcmp(sigSesAndUnits(:,1),sigSessions{ses})) , 2);
    
            for u = 1:length(thisSesSigUnits) % for each significant unit in this session
                subPltCount = subPltCount+1;
                % read out unit's spike counts, already grouped
                eval(['output = spkDataArr' num2str(areA) '.(thisSesSigUnits{u});'])
                
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
                
                % get the parameter values for each group
                paraVal = cell2mat({output.Value}');
        
                % find the right place to plot the tuning curve!
                if subPltCount > (numRows*numColumns)
                    saveas(fg,[SpkCountsPath 'Results\Active Units\' Monkey '_' upper(ParaNames{sigPara}) '_SigUnitsTuning_' Area{areA} '_' num2str(fgCount)], 'fig')
                    saveas(fg,[SpkCountsPath 'Results\Active Units\' Monkey '_' upper(ParaNames{sigPara}) '_SigUnitsTuning_' Area{areA} '_' num2str(fgCount)], 'bmp')
                    close(fg)
                    clear fg
                    fg = figure('Name',[Monkey '_' upper(ParaNames{sigPara}) ' Significant Units_' alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']'], 'Position', [40 70 0 0] + 0.8*get(0,'ScreenSize'));
                    fgCount = fgCount+1;
                    subPltCount = 1;
                end            
                figure(fg), 
                subplot(numRows, numColumns, subPltCount), hold on, box on
                xlabel('EV', 'FontWeight','Bold'); 
                ylabel('Mean Firing Rate', 'FontWeight','Bold')
                title([ {[alignedEvent(2:end-1) '[' num2str(timeInterval(1)) ',' num2str(timeInterval(2)) ']']},...
                        regexprep(sigSessions(ses),'_','-') ,...
                        regexprep(thisSesSigUnits(u),'_','-') ], ...
                        'FontWeight','Normal', 'FontSize',8)
       
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
        
                % add GLM coefficients for this unit
                eval(['aa = find( strcmp(allUnitIDsArray' num2str(areA) '(:,1), sigSessions(ses)) .* strcmp(allUnitIDsArray' num2str(areA) '(:,2), thisSesSigUnits(u)))'])
                eval(['bGLM = bGLMUnitArr' num2str(areA) ';'])
                switch ParaNames{sigPara}
                    case 'ev'
                        txtGLM = [{'GLM coeff.'} , {'-------'} , {['\color[rgb]{1 0 0}bEV = ' num2str(bGLM(aa,2))]} , {['\color[rgb]{0 0 0}bVar = ' num2str(bGLM(aa,3))]} , {['\color[rgb]{0 0 0}bEV*Var = ' num2str(bGLM(aa,4))]}];
                    case 'var'
                        txtGLM = [{'GLM coeff.'} , {'-------'} , {['\color[rgb]{0 0 0}bEV = ' num2str(bGLM(aa,2))]} , {['\color[rgb]{1 0 0}bVar = ' num2str(bGLM(aa,3))]} , {['\color[rgb]{0 0 0}bEV*Var = ' num2str(bGLM(aa,4))]}];
                    case 'int'
                        txtGLM = [{'GLM coeff.'} , {'-------'} , {['\color[rgb]{0 0 0}bEV = ' num2str(bGLM(aa,2))]} , {['\color[rgb]{0 0 0}bVar = ' num2str(bGLM(aa,3))]} , {['\color[rgb]{1 0 0}bEV*Var = ' num2str(bGLM(aa,4))]}];
                end
                text(2.7,min(ylim)+(0.2*diff(ylim)),txtGLM)
                clear aa txtGLM
        
                spkCnt = cell2mat({output.SpikeCount}');
                paraVal = cell2mat(  cellfun(@(x,y)  repmat(x, length(y), 1), {output.Value}',  {output.SpikeCount}', 'UniformOutput', 0)  );
        
                [p,~,~] = anovan(spkCnt, paraVal, 'model', 'interaction', 'varnames', {'EV','Var'}, 'display', 'off');
                txtANOVA = [{'ANOVA'} , {'---------'} , {['pEV = ' num2str(p(1))]} , {['pVar = ' num2str(p(2))]} , {['pEV*Var = ' num2str(p(3))]}];
                text(0.3,min(ylim)+(0.2*diff(ylim)),txtANOVA)
                clear txtANOVA p spkCnt paraVal output
        
            end % end of thisSesSigUnits
    
        end % end of sigSessions
    
        saveas( fg, [SpkCountsPath 'Results\Active Units\' Monkey '_' upper(ParaNames{sigPara}) '_SigUnitsTuning_' Area{areA} '_' num2str(fgCount)], 'fig' )
        saveas( fg, [SpkCountsPath 'Results\Active Units\' Monkey '_' upper(ParaNames{sigPara}) '_SigUnitsTuning_' Area{areA} '_' num2str(fgCount)], 'bmp' )
        close(fg)
        clear fg
    end % end of sigPara<...

    clear allUnit* pGLM* bGLM* ev* var* int* cue*

end % end of areA
    
    
    