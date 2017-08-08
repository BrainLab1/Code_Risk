

%% zscore reaction times for each session; collect the grouped RTs and lick signal from all sessions 

allSesZscrSuccessRT = cell(length(allSesGroupedData{1}), 1);
allSesSuccessLick = cell(length(allSesGroupedData{1}), 1);

for ses = 1:length(allSesGroupedData)
    aa = allSesGroupedData{ses};
    numSuccessInGrp = [];
    sesSuccessRT = []; % size(.)= numSuccessTrialsInSession x 3-> column1:reaction time of the trial; column2:expectedReward; column3:Variance; 
    sesSuccessLick = []; % size(.)= numSuccessTrialsInSession x 3-> column1:reaction time of the trial; column2:expectedReward; column3:Variance; 
    expReward = [];
    rewardVar = [];
    for gr=1:length(aa)
        numSuccessInGrp(gr,1) = length(aa(gr).SuccessRT);
        expReward(gr) = str2num(aa(gr).Value(2));
        switch aa(gr).Value(end)
            case 'H'
                rewardVar(gr) = 4;
            case 'L'
                rewardVar(gr) = 1;
            otherwise
                rewardVar(gr) = 0;
        end
        sesSuccessRT = [ sesSuccessRT; [aa(gr).SuccessRT, ...
                                        repmat(expReward(gr), length(aa(gr).SuccessRT), 1), ...
                                        repmat(rewardVar(gr), length(aa(gr).SuccessRT), 1)] ];
        sesSuccessLick = [ sesSuccessLick; [aa(gr).SuccessLick, ...
                                            repmat(expReward(gr), length(aa(gr).SuccessLick), 1), ...
                                            repmat(rewardVar(gr), length(aa(gr).SuccessLick), 1)] ];
    end
    bb = mat2cell([zscore(sesSuccessRT(:,1)), sesSuccessRT(:,2:3)], numSuccessInGrp', 3);
    cc = mat2cell([sesSuccessLick(:,1), sesSuccessRT(:,2:3)], numSuccessInGrp', 3);
    for gr=1:length(aa)
        aa(gr).ZScrSuccessRT = bb{gr};
        allSesZscrSuccessRT(gr)   = {[allSesZscrSuccessRT{gr};   aa(gr).ZScrSuccessRT]};
        allSesSuccessLick(gr)   = {[allSesSuccessLick{gr};   cc{gr}]};
    end
    allSesGroupedData{ses} = aa;
    clear aa bb cc
end
clear ses sesSuccessRT sesSuccessLick

%% remove reaction time outliers
allSesZscrSuccessRT = cell2mat(allSesZscrSuccessRT); % [reactionTime expevtedReward rewardVariance]
numOriginalTrials = size(allSesZscrSuccessRT,1);
% find trials with RT out out of the predefined boundry
trlIdx = find( (allSesZscrSuccessRT(:,1) < -rtOutlierFactor) + (allSesZscrSuccessRT(:,1) > rtOutlierFactor) );
allSesZscrSuccessRT(trlIdx,:) = [];
numLostTrialsRT = length(trlIdx);
clear trlIdx

%% plot figure 5: reaction times
[~,tblRT] = anovan(allSesZscrSuccessRT(:,1), {allSesZscrSuccessRT(:,2), allSesZscrSuccessRT(:,3)}, 'model', 'interaction', 'varnames',{'ExpectedValue','Variance'});
% find trials with [ExpReward Variance] 
expReward = [3 6 9];
rewVar = [0 1 4];
for vr=1:3
    trlIdx = [];
    for ex = 1:3
        trlIdx = find( (allSesZscrSuccessRT(:,2) == expReward(ex)) .* (allSesZscrSuccessRT(:,3)==rewVar(vr)) );
        eval(['crv' num2str(rewVar(vr)) '(ex) = mean(allSesZscrSuccessRT(trlIdx,1));'])
        eval(['err' num2str(rewVar(vr)) '(ex) = std(allSesZscrSuccessRT(trlIdx,1))/sqrt(length(trlIdx));'])
        clear trlIdx
    end
end
clear vr ex

fg5 = figure('Name', monkeyName);
subplot(2,1,1), hold on, box on
title(['Outlier criterion: ' num2str(rtOutlierFactor) '*std.   ' num2str(numLostTrialsRT) ' out of ' num2str(numOriginalTrials) ' trials removed'])
xlabel('Expected value')
ylabel({'Ave. reaction time','+/- SEM'})
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

%% plot figure 5: licking probability
allSesSuccessLick = cell2mat(allSesSuccessLick);  % [lickingSignal expevtedReward rewardVariance]

[~,tblLick] = anovan(allSesSuccessLick(:,1), {allSesSuccessLick(:,2), allSesSuccessLick(:,3)}, 'model', 'interaction', 'varnames',{'ExpectedValue','Variance'});
% find trials with [ExpReward Variance] 
expReward = [3 6 9];
rewVar = [0 1 4];
for vr=1:3
    trlIdx = [];
    for ex = 1:3
        trlIdx = find( (allSesSuccessLick(:,2) == expReward(ex)) .* (allSesSuccessLick(:,3)==rewVar(vr)) );
        eval(['crv' num2str(rewVar(vr)) '(ex) = mean(allSesSuccessLick(trlIdx,1));'])
        eval(['err' num2str(rewVar(vr)) '(ex) = std(allSesSuccessLick(trlIdx,1))/sqrt(length(trlIdx));'])
    end
end
clear vr ex

figure(fg5)
subplot(2,1,2), hold on, box on
title(['Mean lick signal at (' lickTimeRelativeEvent ' [' num2str(lickBestTime(1)) ' ' num2str(lickBestTime(2)) '])'])
xlabel('Expected value')
ylabel({'Ave. licking signal','+/- SEM'})
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

