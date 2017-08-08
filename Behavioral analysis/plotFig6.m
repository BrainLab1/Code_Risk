
%% zscore reaction times for each session; collect the grouped RTs and lick signal from all sessions 

allSesSucRTVar = cell(length(allSesGroupedData{1}), 1);
allSesSucLickVar = cell(length(allSesGroupedData{1}), 1);

for ses = 1:length(allSesGroupedData)
    aa = allSesGroupedData{ses};
    numSuccessInGrp = [];
    sesSuccessRT = []; % size(.)= numSuccessTrialsInSession x 3-> column1:reaction time of the trial; column2:expectedReward; column3:Variance; 
    sesSucRTVar = []; % size(.)= numSuccessTrialsInSession x 3-> column1:reaction time of the trial; column2:expectedReward; column3:Variance; 
    sesSucLickVar = []; % size(.)= numSuccessTrialsInSession x 3-> column1:reaction time of the trial; column2:expectedReward; column3:Variance; 
    expReward = [];
    rewardVar = [];
    for gr=1:length(aa)
        numSuccessInGrp(gr,1) = length(aa(gr).SuccessRT);
        expReward(gr,1) = str2num(aa(gr).Value(2));
        switch aa(gr).Value(end)
            case 'H'
                rewardVar(gr,1) = 4;
            case 'L'
                rewardVar(gr,1) = 1;
            otherwise
                rewardVar(gr,1) = 0;
        end
        sesSuccessRT = [ sesSuccessRT; [aa(gr).SuccessRT, ...
                                        repmat(expReward(gr), length(aa(gr).SuccessRT), 1), ...
                                        repmat(rewardVar(gr), length(aa(gr).SuccessRT), 1)] ];
        sesSucLickVar = [ sesSucLickVar; [var(aa(gr).SuccessLick), expReward(gr), rewardVar(gr)] ];
    end
    clear gr
    % zscore the session reaction times
    bb = mat2cell([zscore(sesSuccessRT(:,1)), sesSuccessRT(:,2:3)], numSuccessInGrp, 3);
    for gr=1:length(aa)
        % remove RT outliers: for each group, find trials with RT out out of the predefined boundry
        tmpbb =  bb{gr};
        trlIdx = find( (tmpbb(:,1) < -rtOutlierFactor) + (tmpbb(:,1) > rtOutlierFactor) );
        tmpbb(trlIdx,:) = [];
        clear trlIdx
        % add the group zscored reaction time and lick to the population group data
        aa(gr).ZScrSuccessRT = tmpbb;
        sesSucRTVar(gr,:)      = [var(aa(gr).ZScrSuccessRT(:,1)), unique(tmpbb(:,2)), unique(tmpbb(:,3))];
        aa(gr).GroupZScrSucRTVar  = sesSucRTVar(gr,:);
        aa(gr).GroupSucLickVar    = sesSucLickVar(gr,:);
        allSesSucRTVar(gr)     = {[allSesSucRTVar{gr};     sesSucRTVar(gr,:)]};
        allSesSucLickVar(gr)   = {[allSesSucLickVar{gr};   sesSucLickVar(gr,:)]};
    end
    allSesGroupedData{ses} = aa;
    clear aa bb cc
end
clear ses sesSuccessRT sesSuccessLick

%% plot figure 6: reaction times
allSesSucRTVar = cell2mat(allSesSucRTVar);

[~,tblRT] = anovan(allSesSucRTVar(:,1), {allSesSucRTVar(:,2), allSesSucRTVar(:,3)}, 'model', 'interaction', 'varnames',{'ExpectedValue','Variance'});
% find sessions with [ExpReward Variance] 
expReward = [3 6 9];
rewVar = [0 1 4];
for vr=1:3
    sesIdx = [];
    for ex = 1:3
        sesIdx = find( (allSesSucRTVar(:,2) == expReward(ex)) .* (allSesSucRTVar(:,3)==rewVar(vr)) );
        eval(['crv' num2str(rewVar(vr)) '(ex) = mean(allSesSucRTVar(sesIdx,1));'])
        eval(['err' num2str(rewVar(vr)) '(ex) = std(allSesSucRTVar(sesIdx,1))/sqrt(length(sesIdx));'])
        clear sesIdx
    end
end
clear vr ex

fg6 = figure('Name', monkeyName);
subplot(2,1,1), hold on, box on
title(['Outlier criterion: ' num2str(rtOutlierFactor) '*std'])
xlabel('Expected value')
ylabel({'Ave. reaction time variance','+/- SEM'})
txt = {};
for vr=1:3
    eval(['errorbar( 1:3, crv' num2str(rewVar(vr)) ', err' num2str(rewVar(vr)) ', ''Color'', colorCode(vr+3,:), ''LineWidth'', 2.5)'])
    txt = [txt; {['\color[rgb]{' num2str(colorCode(vr+3,:)) '} variance = ' num2str(rewVar(vr)) ]}];
end
txt = [txt; {'\color[rgb]{0 0 0}  '}];
txt = [txt; {['p' tblRT{2,1} ' = ' num2str(tblRT{2,7})]}];
txt = [txt; {['p' tblRT{3,1} ' = ' num2str(tblRT{3,7})]}];
txt = [txt; {['pInteraction  = ' num2str(tblRT{4,7})]}];

set(gca, 'XLim', [0 4], 'XTick', [1 2 3], 'XTickLabel', {'3','6','9'})
text(mean(xlim), mean(ylim), txt)

clear crv0 err0 crv1 err1 crv4 err4 txt

%% plot figure 6: licking probability
allSesSucLickVar = cell2mat(allSesSucLickVar);  % [lickingSignal expevtedReward rewardVariance]

[~,tblLick] = anovan(allSesSucLickVar(:,1), {allSesSucLickVar(:,2), allSesSucLickVar(:,3)}, 'model', 'interaction', 'varnames',{'ExpectedValue','Variance'});
% find trials with [ExpReward Variance] 
expReward = [3 6 9];
rewVar = [0 1 4];
for vr=1:3
    trlIdx = [];
    for ex = 1:3
        sesIdx = find( (allSesSucLickVar(:,2) == expReward(ex)) .* (allSesSucLickVar(:,3)==rewVar(vr)) );
        eval(['crv' num2str(rewVar(vr)) '(ex) = mean(allSesSucLickVar(sesIdx,1));'])
        eval(['err' num2str(rewVar(vr)) '(ex) = std(allSesSucLickVar(sesIdx,1))/sqrt(length(sesIdx));'])
        clear sesIdx
    end
end
clear vr ex

figure(fg6)
subplot(2,1,2), hold on, box on
title(['Licking signal variance at (' lickTimeRelativeEvent ' [' num2str(lickBestTime(1)) ' ' num2str(lickBestTime(2)) '])'])
xlabel('Expected value')
ylabel({'Ave. licking signal variance','+/- SEM'})
txt = {};
for vr=1:3
    eval(['errorbar( 1:3, crv' num2str(rewVar(vr)) ', err' num2str(rewVar(vr)) ', ''Color'', colorCode(vr+3,:), ''LineWidth'', 2.5)'])
    txt = [txt; {['\color[rgb]{' num2str(colorCode(vr+3,:)) '} variance = ' num2str(rewVar(vr)) ]}];
end
txt = [txt; {'\color[rgb]{0 0 0}  '}];
txt = [txt; {['p' tblLick{2,1} ' = ' num2str(tblLick{2,7})]}];
txt = [txt; {['p' tblLick{3,1} ' = ' num2str(tblLick{3,7})]}];
txt = [txt; {['pInteraction  = ' num2str(tblLick{4,7})]}];

set(gca, 'XLim', [0 4], 'XTick', [1 2 3], 'XTickLabel', {'3','6','9'})
text(mean(xlim), mean(ylim), txt)

