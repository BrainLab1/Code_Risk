% 08.05.2018: this is the same scripts as R13_ApplyDPCA but with the new channel-electrode map 
% 10.03.2018

allUnitsPSTH_Arr_Base{1,2} = [1 1]; %[3,0];
allUnitsPSTH_Arr_Base{2,2} = [1 2]; %[3,1];
allUnitsPSTH_Arr_Base{3,2} = [1 3]; %[3,4];
allUnitsPSTH_Arr_Base{4,2} = [2 1]; %[6,0];
allUnitsPSTH_Arr_Base{5,2} = [2 2]; %[6,1];
allUnitsPSTH_Arr_Base{6,2} = [2 3]; %[6,4];
allUnitsPSTH_Arr_Base{7,2} = [3 1]; %[9,0];
allUnitsPSTH_Arr_Base{8,2} = [3 2]; %[9,1];
allUnitsPSTH_Arr_Base{9,2} = [3 3]; %[9,4];

aa = [];
for un=1:size(allUnitsPSTH_Arr_Base{1,1},1)
    for c = 1:size(allUnitsPSTH_Arr_Base,1)
        aa(un, allUnitsPSTH_Arr_Base{c,2}(1), allUnitsPSTH_Arr_Base{c,2}(2), :) = allUnitsPSTH_Arr_Base{c,1}(un,:);
    end
end
clear c un

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'EV', 'Var', 'Condition-independent', 'EV/Var Interaction'};

[W,V,whichMarg] = dpca(aa, 20, 'combinedParams', combinedParams);

% centering
X = aa(:,:);
aa_mean = bsxfun(@minus, aa, mean(X,2));
X = bsxfun(@minus, X, mean(X,2));

% % % Z = W'*X;
Z = {};
for ev = 1:size(aa,2)
    for vr = 1:size(aa,3)
        Z{ev,vr} = W' * squeeze(aa(:,ev,vr,:));
    end
end
clear ev vr

explVar = dpca_explainedVariance(aa, W, V, 'combinedParams', combinedParams);

timeEvents = [0 0.4 0.8]; % event times in sec to be shown on the figures
margLabelColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

dpca_plot(aa, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margLabelColours, ...
    'whichMarg', whichMarg,                 ...
    'time', newT,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);

%% plot dPCA
if 0 
pcs2plot = [1:4];
numRows = 2; 
numColumns = 2;
figure('Name', [Monkey ', ' Area])
for pcs = 1:length(pcs2plot)
    subplot(numRows, numColumns, pcs), hold on, box on
    title(['PC ' num2str(pcs2plot(pcs))])
    ccnt = 0;
    for ev = 1:3
        for vr = 1:3
           ccnt = ccnt+1;
           plot(newT, Z{ev,vr}(pcs2plot(pcs),:), 'Color', colorCode(ccnt,:))
        end
    end
    legend('[3 0]','[3 1]','[3 4]','[6 0]','[6 1]','[6 4]','[9 0]','[9 1]','[9 4]')
end

end