

clear all;
clc;
dbstop if error

% Set the parameters bellow so that [main_folder data_folder bhv_folder] 
% points to the directory where all the preprocessed files are stored.
main_folder          = 'Z:\'; 
data_folder          = 'data\'; % 
original_data_folder = 'Risk\Original Data_Extracted\';
save_dir             = 'X:\Report 2\Bahareh\BehaviorData\';

colorCode = [  204 204 255  % F3
               255 153 153  % F6
               153 255 204   % F9
               102 102 255   % B3L
               255 51  51    % B6L
               51  255 102   % B9L              
               0   0   255   % B3H
               204 0   0     % B6H
               0  153 153]/255; % B9H   
           
rtOutlierFactor = 2; % outliers are trials out of this intervel [-1 1]*trOutlierFactor*Std; Std for each monkey is different
monkeyName = 'MacDuff';  % 'Moj' or 'Mac'
groupingCriteria       = 'type & Congruence' ; % this parameter shows based on what criteria trials are grouped.
Fs = 1000;   % assume 1K sampling rate

%% Read out list of all the files related to this session
dataPath = [main_folder data_folder original_data_folder];
allFiles = dir(dataPath);    

% remove the '.' and '..' from the folders list and filter sessions for monkey names
idx = [];
for ses = 1:numel(allFiles) % for each session
    if (numel(allFiles(ses).name) < 3) || ~(strcmp(allFiles(ses).name(1:3), monkeyName(1:3)))  
        idx = [idx; ses];
    end
end
allFiles(idx) = [];
clear ses idx

%% Extract and group licking signal and rea
allSesGroupedData = {};

for ses = 1:numel(allFiles) % for each session
    cfg = [];
    sessionFolder = allFiles(ses).name;
 if 0 % ---------------------------------------------------------------------------------------------------------------------   
    bhvFileInfo = dir([dataPath sessionFolder '\*.bhv']); 
    bhvFilePath = [dataPath sessionFolder '\' bhvFileInfo.name];
        
    lfpFileInfo    = dir([dataPath sessionFolder '\*.ns2']);   
    if length(lfpFileInfo) > 1
        warning(['Multiple NS files in folder ' [dataPath sessionFolder]]);
        continue;
    end       
    lfpFilePath = [dataPath sessionFolder '\' lfpFileInfo.name];

    % form the configuration structure for 'ft_definetrial()'
    tmp_cfg                 = [];
    tmp_cfg.headerformat    = 'blackrock_nsx';
    tmp_cfg.data_dir        = dataPath;
    tmp_cfg.trialfun        = 'fieldtrip_trialfun_RiskBhv'; % string with function name, see below (default = 'ft_trialfun_general')
    tmp_cfg.trialdef.pre    = 1; % time in seconds; this only indicates absolut value, do not use '-'
    tmp_cfg.trialdef.post   = 1; % time in seconds
        
    tmp_cfg.headerfile  = lfpFilePath;
    tmp_cfg.dataset     = lfpFilePath;
    tmp_cfg.session_dir = sessionFolder;
        
    % Config for aligning to cue
    tmp_cfg.trialdef.eventtype  = alignEvent;
    tmp_cfg.trialdef.interval   = evetTimeInterval;

    % read out the data the same way it is done in ReadRiskDataBhv.m
    new_cfg = ft_definetrial(tmp_cfg);
    
    % save the extracted data on the hard drive for future use
    save([save_dir allFiles(ses).name '(' alignEvent ')'], 'new_cfg', 'alignEvent', 'evetTimeInterval')
 end % ----------------------------------------------------------------------------------------------------------------------

    load ([save_dir allFiles(ses).name '(cue)'])
 
    % group trials for 9 task conditions
    [output] = GroupTrials(new_cfg.event, groupingCriteria);

    % convert the event tructure to table
    eventTable = struct2table(new_cfg.event);

    % get the trial index for successful trials per condition
    for gr = 1:length(output)
        output(gr).SuccessTrlIdx = output(gr).TrialIdx( eventTable.TrialErrorCode(output(gr).TrialIdx) == 0 );
        % for successful trials of this group, collect the RT 
        tempRT = [];
        for tr = 1:length(output(gr).SuccessTrlIdx)
            trIdx = output(gr).SuccessTrlIdx(tr);
            if (ses==9 && trIdx==133)
                continue
            end
            % read out the reaction time for this successful trial of this group
            tempRT = [tempRT; eventTable.ReactionTime(trIdx)];
            clear trIdx
        end
        output(gr).SuccessRT   = tempRT; 
        clear  tempRT tr 

    end
    
    allSesGroupedData{ses} = output;
    clear gr output eventTable new_cfg tmp_cfg temp 
end
clear ses

%% zscore reaction times of successful trials, for each session independently
allSesZscSucRT = cell(length(allSesGroupedData{1}), 1); % cell array with the same size as the num. of groups; in each cell all trials of the group from all sessions are included, plus expVal and Variance  
for ses = 1:numel(allSesGroupedData) % for each session
    sesSucRT = [];
    numTrlPerGrp = [];
    aa = allSesGroupedData{ses};
    for gr = 1:length(allSesGroupedData{ses})
        sesSucRT = [sesSucRT; aa(gr).SuccessRT];
        numTrlPerGrp = [numTrlPerGrp; length(aa(gr).SuccessRT)];
    end
    clear gr
    % zscore all the successful RTs for this session
    sesZscSucRT = zscore(sesSucRT);
    % convert re-group the RTs
    sesZscSucRT = mat2cell(sesZscSucRT, numTrlPerGrp, 1);
    % add a new field to this session data for zscored RTs
    for gr = 1:length(allSesGroupedData{ses})
        aa(gr).ZscrSuccessRT = sesZscSucRT{gr};
        expVal  = str2num(aa(gr).Value(2));
        congIdx = str2num(aa(gr).Value(end-1:end));
        switch aa(gr).Value(3)
            case 'H'
                rewVar = 4;
            case 'L'
                rewVar = 1;
            otherwise
                rewVar = 0;
        end
        allSesZscSucRT{gr} = [ allSesZscSucRT{gr}; [ sesZscSucRT{gr} repmat([expVal rewVar congIdx ses], length(sesZscSucRT{gr}), 1) ] ];
        clear expVal congIdx rewVar
    end
    clear gr
    allSesGroupedData{ses} = aa;
    clear sesZscSucRT numTrlPerGrp sesSucRT aa
end
clear ses

%% remove reaction time outliers
allSesZscSucRT = cell2mat(allSesZscSucRT); % [reactionTime expevtedReward rewardVariance congruencyIdx sessionIdx]
numOriginalTrials = size(allSesZscSucRT,1);
% find trials with RT out out of the predefined boundry
trlIdx = find( (allSesZscSucRT(:,1) < -rtOutlierFactor) + (allSesZscSucRT(:,1) > rtOutlierFactor) );
allSesZscSucRT(trlIdx,:) = [];
numLostTrialsRT = length(trlIdx);
clear trlIdx

%%  plot fig7_1
[~,tblRT] = anovan(allSesZscSucRT(:,1), {allSesZscSucRT(:,2), allSesZscSucRT(:,3), allSesZscSucRT(:,4)}, 'model', 'interaction', 'varnames',{'ExpectedValue','Variance','Congruence'});
% find trials with [ExpReward Variance] 
expReward = [3 6 9];
rewVar = [0 1 4];
congIdx = [-1 1];
for vr=1:3
    trlIdx = [];
    for ex = 1:3
        % for this variance and expValue, find congruent trials
        trlIdx = find( (allSesZscSucRT(:,2) == expReward(ex)) .* (allSesZscSucRT(:,3)==rewVar(vr)) .* (allSesZscSucRT(:,4)==1) );
        eval(['crv' num2str(rewVar(vr)) 'Cong(ex) = mean(allSesZscSucRT(trlIdx,1));'])
        eval(['err' num2str(rewVar(vr)) 'Cong(ex) = std(allSesZscSucRT(trlIdx,1))/sqrt(length(trlIdx));'])

        trlIdx = [];
        % for this variance and expValue, find incongruent trials
        trlIdx = find( (allSesZscSucRT(:,2) == expReward(ex)) .* (allSesZscSucRT(:,3)==rewVar(vr)) .* (allSesZscSucRT(:,4)==-1) );
        eval(['crv' num2str(rewVar(vr)) 'Incong(ex) = mean(allSesZscSucRT(trlIdx,1));'])
        eval(['err' num2str(rewVar(vr)) 'Incong(ex) = std(allSesZscSucRT(trlIdx,1))/sqrt(length(trlIdx));'])
        clear trlIdx
    end
end
clear vr ex

fg = figure('Name', monkeyName);
subplot(2,1,1), hold on, box on
title(['Outlier criterion: ' num2str(rtOutlierFactor) '*std.   ' num2str(numLostTrialsRT) ' out of ' num2str(numOriginalTrials) ' trials removed'])
xlabel('Expected value')
ylabel({'Ave. reaction time','+/- SEM'})
txt = {};
for vr=1:3
    eval(['errorbar( 1:3, crv' num2str(rewVar(vr)) 'Cong,   err' num2str(rewVar(vr)) 'Cong,   ''Color'', colorCode(vr+3,:), ''LineWidth'', 2.5)'])
    eval(['errorbar( 1:3, crv' num2str(rewVar(vr)) 'Incong, err' num2str(rewVar(vr)) 'Incong, ''Color'', colorCode(vr+3,:), ''LineWidth'', 2.5, ''LineStyle'',''--'')'])
    txt = [txt; {['\color[rgb]{' num2str(colorCode(vr+3,:)) '} Variance = ' num2str(rewVar(vr)) ]}];
end
txt = [txt; {'\color[rgb]{0 0 0}  '}];
txt = [txt; {['p' tblRT{2,1} ' = ' num2str(tblRT{2,7})]}];
txt = [txt; {['p' tblRT{3,1} ' = ' num2str(tblRT{3,7})]}];
txt = [txt; {['p' tblRT{4,1} ' = ' num2str(tblRT{4,7})]}];
txt = [txt; {['p ' tblRT{5,1} ' = ' num2str(tblRT{5,7})]}];
txt = [txt; {['p ' tblRT{6,1} ' = ' num2str(tblRT{6,7})]}];
txt = [txt; {['p ' tblRT{7,1} ' = ' num2str(tblRT{7,7})]}];

set(gca, 'XLim', [0 4], 'XTick', [1 2 3], 'XTickLabel', {'3','6','9'})
text(mean(xlim), mean(ylim), txt)
clear crv* err* txt

%%  plot fig7_2
numSession = length(unique(allSesZscSucRT(:,5)));
expReward = [3 6 9];
rewVar = [0 1 4];
congIdx = [-1 1];

for ses = 1:numSession
    for vr=1:3
        trlIdx = [];
        for ex = 1:3
            % for this variance and expValue, find congruent trials for this session
            trlIdx = find( (allSesZscSucRT(:,2) == expReward(ex)) .* (allSesZscSucRT(:,3)==rewVar(vr)) .* (allSesZscSucRT(:,4)==1) .* (allSesZscSucRT(:,5)==ses) );
            eval(['crv' num2str(rewVar(vr)) 'Cong(ses,ex) = mean(allSesZscSucRT(trlIdx,1));'])

            trlIdx = [];
            % for this variance and expValue, find incongruent trials
            trlIdx = find( (allSesZscSucRT(:,2) == expReward(ex)) .* (allSesZscSucRT(:,3)==rewVar(vr)) .* (allSesZscSucRT(:,4)==-1) .* (allSesZscSucRT(:,5)==ses) );
            eval(['crv' num2str(rewVar(vr)) 'Incong(ses,ex) = mean(allSesZscSucRT(trlIdx,1));'])
            clear trlIdx
        end
        eval(['crv' num2str(rewVar(vr)) 'Delta(ses,:) =  crv' num2str(rewVar(vr)) 'Cong(ses,:) - crv' num2str(rewVar(vr)) 'Incong(ses,:) ;'])        
    end
end
clear vr ex

figure(fg);
subplot(2,1,2), hold on, box on
xlabel('Expected value')
ylabel({'congRT - incongRT','+/- SEM'})
txt = {};
for vr=1:3
    eval(['err = std(crv' num2str(rewVar(vr)) 'Delta)./sqrt(size(crv' num2str(rewVar(vr)) 'Delta,1))'])
    eval(['errorbar( 1:3, mean(crv' num2str(rewVar(vr)) 'Delta),   err,   ''Color'', colorCode(vr+3,:), ''LineWidth'', 2.5)'])
    txt = [txt; {['\color[rgb]{' num2str(colorCode(vr+3,:)) '} Variance = ' num2str(rewVar(vr)) ]}];
end
% % txt = [txt; {'\color[rgb]{0 0 0}  '}];
% % txt = [txt; {['p' tblRT{2,1} ' = ' num2str(tblRT{2,7})]}];
% % txt = [txt; {['p' tblRT{3,1} ' = ' num2str(tblRT{3,7})]}];
% % txt = [txt; {['p' tblRT{4,1} ' = ' num2str(tblRT{4,7})]}];
% % txt = [txt; {['p ' tblRT{5,1} ' = ' num2str(tblRT{5,7})]}];
% % txt = [txt; {['p ' tblRT{6,1} ' = ' num2str(tblRT{6,7})]}];
% % txt = [txt; {['p ' tblRT{7,1} ' = ' num2str(tblRT{7,7})]}];

set(gca, 'XLim', [0 4], 'XTick', [1 2 3], 'XTickLabel', {'3','6','9'})
text(mean(xlim), mean(ylim), txt)
clear crv* err* txt

