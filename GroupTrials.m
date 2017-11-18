% last update: 18.11.2017 by Saeed: new conditions for uncertainty task is added
% last update: 17.09.2017 by Bahareh: new condition added for 'CuePos & TargetPos & CurrEV & CurrVar' 
% last update: 16.09.2017 by Bahareh: new condition added for 'PreTrlVar & PreTrlOutcome & CurrTrlEV' 
% last update: 04.09.2017 by Bahareh: new condition added for 'PreTrlEV & PreTrlOutcome & CurrTrlEV' 
% last update: 29.08.2017 by Bahareh: new condition added for 'PreTrlEV & PreTrlVar & CurrTrlEV' 
% last update: 29.08.2017 by Bahareh: new condition added for 'expecter_reward & RewardVariance & Outcome'
% last update: 28.08.2017 by Bahareh: new condition added for 'expecter_reward & RewardVariance'
% last update: 20.08.2017 by Bahareh: new condition added for 'PreTrlOutcome & CurrentTrlVar' where successful trials for which the 
%      previos trial was rewarded, are grouped for reward variance of the current trial and outcome (win or lose) of the previouse trial
% last update: 20.08.2017 by Bahareh: new condition added for 'PreTrlOutcome & CurrentTrlEV' where successful trials for which the 
%      previos trial was rewarded, are grouped for expected value of the current trial and outcome (win or lose) of the previouse trial
% last update: 19.08.2017 by Saeed: new condition added for 'OutlierIndex'
% last update: 19.08.2017 by Bahareh: new condition added for 'expected_reward & Congruence'
% last update: 15.08.2017 by Saeed: new conditions added for 'PreTrialEV', 'PreTrialVAR' & 'PreTrialEV&VAR'
% last update: 30.07.2017 by Bahareh: new conditions added for 'type & Congruence', 
% last update: 25.07.2017 by Bahareh: new conditions added for '18 task conditions', 'SaccadeLaterality', 
% last update: 24.07.2017 by Bahareh: new conditions added for 'CueTargetLaterality' and 'CueTargetCongruency' and 'Entropy'
% last update by Bahareh 22.07.2017: a new condition was added for 'Variance'
% created by Bahareh on 18.07.2017
% inputs: 
%   - event: is a cell array of structures. it is created by fieldtrip_trialfun_Risk.m or fieldtrip_trialfun_RiskBhv.m    
%   - grType: indicates how to group the trials. it could be any of the fields of the event structures 
% 
% 
% output:
%   - output: cell array


function  [output] = GroupTrials(event, grType)

switch grType
    case 'type'
        eventTable = struct2table(event);
        possibleVals = unique(eventTable.(grType));
        output = [];
        for i=1:length(possibleVals)
            output = [output; struct('TrialIdx', find(strcmp(eventTable.(grType), possibleVals{i})), ...
                                     'Value', possibleVals{i}, ...
                                     'GroupingType', grType)];
        end
        
    case {'value', 'expected_reward', 'cue_pos', 'target_pos', 'pre_good', 'taskCondCode','subjectID','SessionID','OutlierIndex'}
        eventTable = struct2table(event);
        possibleVals = unique(eventTable.(grType));
        output = [];
        for i=1:length(possibleVals)
            output = [output; struct('TrialIdx', find(eventTable.(grType) == possibleVals(i)), ...
                                     'Value', possibleVals(i), ...
                                     'GroupingType', grType)];
        end
        
    case 'Variance'
        eventTable = struct2table(event);
        output = [];
        % get the trial indices for variance 0
        idx = [];
        for tr = 1:length(eventTable.type)
            if strcmp(eventTable.type{tr}(1) , 'F')
                idx = [idx; tr];
            end
        end
        clear tr
        output = [output; struct('TrialIdx', idx, ...
                                 'Value', 0, ...
                                 'GroupingType', grType)];
        clear idx
        
        % get the trial indices for variance 4
        idx = [];
        for tr = 1:length(eventTable.type)
            if ( strcmp(eventTable.type{tr}(1) , 'B') && strcmp(eventTable.type{tr}(end) , 'L') )
                idx = [idx; tr];
            end
        end
        clear tr
        output = [output; struct('TrialIdx', idx, ...
                                 'Value', 1, ...
                                 'GroupingType', grType)];
        clear idx

        % get the trial indices for variance 4
        idx = [];
        for tr = 1:length(eventTable.type)
            if ( strcmp(eventTable.type{tr}(1) , 'B') && strcmp(eventTable.type{tr}(end) , 'H') )
                idx = [idx; tr];
            end
        end
        clear tr
        output = [output; struct('TrialIdx', idx, ...
                                 'Value', 4, ...
                                 'GroupingType', grType)];
        clear idx
        
    case 'Risk'
        eventTable = struct2table(event);
        allType = {'F','B3L',{'B3H','B6L'},'B9L','B6H','B9H'};
        allRisk = [0     3         6         9    12    18];
        output = [];
        output = [output; struct('TrialIdx', find( strcmp(eventTable.type, {'F3'}) + strcmp(eventTable.type, {'F6'}) + strcmp(eventTable.type, {'F9'}) ), ...
                                 'Value', allRisk(1), ...
                                 'GroupingType', grType)];
        idx1 = find( strcmp(eventTable.type, {'F3'}) + strcmp(eventTable.type, {'F6'}) + strcmp(eventTable.type, {'F9'}) );                     
        output = [output; struct('TrialIdx', find( strcmp(eventTable.type, {'B3L'}) ), ...
                                 'Value', allRisk(2), ...
                                 'GroupingType', grType)];
        idx2= find( strcmp(eventTable.type, {'B3L'}) );                            
        output = [output; struct('TrialIdx', find( strcmp(eventTable.type, {'B3H'}) + strcmp(eventTable.type, {'B6L'}) ), ...
                                 'Value', allRisk(3), ...
                                 'GroupingType', grType)];
        idx3= find( strcmp(eventTable.type, {'B3H'}) + strcmp(eventTable.type, {'B6L'}) );                     
       for i=4:length(allRisk)
            output = [output; struct('TrialIdx', find( strcmp(eventTable.type , allType(i))), ...
                                     'Value', allRisk(i), ...
                                     'GroupingType', grType)];
            eval(['idx' num2str(i) '=find( strcmp(eventTable.type , allType(i)));'])                     
       end
        
   case 'CueTargetLaterality'
        eventTable = struct2table(event);
        output = [];
        possibleCueVals = [-1 1];
        possibleTarVals = [-1 1];
        positionCode = [1 2; ...
                        3 4]; 
        for cu = 1:length(possibleCueVals)
            for tar = 1:length(possibleTarVals)
                output = [output; struct('TrialIdx', find( (eventTable.cue_pos == possibleCueVals(cu)) .* (eventTable.target_pos == possibleTarVals(tar)) ), ...
                                         'Value', positionCode(cu, tar), ...
                                         'GroupingType', grType, ...
                                         'AdditionalInfo', ['[Cue Target] = [' num2str(possibleCueVals(cu)) ' ' num2str(possibleTarVals(tar)) ']' ])];
            end
        end

   case 'CueTargetCongruency'
        eventTable = struct2table(event);
        output = [];
        output = [output; struct('TrialIdx', find( sign(eventTable.cue_pos .* eventTable.target_pos )>0 ), ...
                                 'Value', 1, ...
                                 'GroupingType', grType)];
        output = [output; struct('TrialIdx', find( sign(eventTable.cue_pos .* eventTable.target_pos )<0 ), ...
                                 'Value', -1, ...
                                 'GroupingType', grType)];

   case 'Entropy'
        output = [];
        idx0 = [];
        idx1 = [];
        for tr = 1:size(event,1)
            switch event(tr).type(1)
                case 'F'
                    idx0 = [idx0; tr];
                case 'B'
                    idx1 = [idx1; tr];
            end
        end
        clear tr
        
        output = [output; struct('TrialIdx', idx0, ...
                                 'Value', 0, ...
                                 'GroupingType', grType)];
        output = [output; struct('TrialIdx', idx1, ...
                                 'Value', 1, ...
                                 'GroupingType', grType)];
                             
   case '18 task conditions' % grop trial by : ExpectedRewardValue & RewardVariance & cue-target Congruency 
        output = [];
        eventTable = struct2table(event);
        congruencyIndx = (eventTable.cue_pos .* eventTable.target_pos) > 0;
        possibleCong   =  unique(congruencyIndx);
        expReward = unique(eventTable.expected_reward);
        rewardVar = unique(eventTable.RewardVariance);
        for ex = 1:length(expReward)
            for vr = 1:length(rewardVar)
                for cg = 1:length(possibleCong)
                    trIdx = find ( (eventTable.expected_reward == expReward(ex)) .* ...
                                   (eventTable.RewardVariance == rewardVar(vr)) .* ...
                                   (congruencyIndx == possibleCong(cg)) );
                
                    output = [output; struct('TrialIdx', trIdx, ...
                                             'Value', [expReward(ex) rewardVar(vr) possibleCong(cg)], ...
                                             'GroupingType', grType)];
                             
                end
            end
        end
        
    case 'SaccadeLaterality'
        eventTable = struct2table(event);
        output = [];
        possibleValues = unique(eventTable.target_pos);
        for i = 1:length(possibleValues)
            output = [output; struct('TrialIdx', find( eventTable.target_pos == possibleValues(i) ), ...
                                     'Value', possibleValues(i), ...
                                     'GroupingType', grType)];
        end
        clear i
        
        
    case 'type & Congruence'
        eventTable = struct2table(event);
        possibleVals = unique(eventTable.type); % this will give you 9 cue types
        possibleVals = [ [possibleVals,  mat2cell(ones(length(possibleVals),1),ones(1,length(possibleVals)),1)];...
                         [possibleVals,  mat2cell(-ones(length(possibleVals),1),ones(1,length(possibleVals)),1)] ];
        output = [];
        congeurncyIdx = eventTable.cue_pos .* eventTable.target_pos;
        for i=1:size(possibleVals,1)
            output = [output; struct('TrialIdx', find( strcmp(eventTable.type, possibleVals{i,1}) .* (congeurncyIdx == possibleVals{i,2}) ), ...
                                     'Value', [possibleVals{i,1} ' ' num2str(possibleVals{i,2})], ...
                                     'GroupingType', grType)];
        end
        
    case 'expected_reward & Congruence'
        eventTable = struct2table(event);
        possibleVals = unique(eventTable.expected_reward); % this will gige you 3 expected reward value
        possibleVals = [ [possibleVals,  ones(length(possibleVals),1)];...
                         [possibleVals, -ones(length(possibleVals),1)] ];
        output = [];
        congeurncyIdx = eventTable.cue_pos .* eventTable.target_pos;
        for i=1:size(possibleVals,1)
            output = [output; struct('TrialIdx', find( (eventTable.expected_reward == possibleVals(i,1)) .* (congeurncyIdx == possibleVals(i,2)) ), ...
                                     'Value', possibleVals(i,:), ...
                                     'GroupingType', grType)];
        end
        
    case 'PreTrlOutcome & CurrentTrlEV'
        output = [];
        win3 = [];  win6 = [];  win9 = [];
        los3 = [];  los6 = [];  los9 = [];
        eventTable = struct2table(event);        
        for tr = 2:size(eventTable,1)  % for each trial
            if ~eventTable.TrialErrorCode(tr) % if current trial was successful
                if ( ~eventTable.TrialErrorCode(tr-1)  && ~strcmp(eventTable.type{tr-1}(1),'F') ) % if previous trial was successful (rewarded) && reawrd variance was not zero   
                    switch eventTable.expected_reward(tr)
                        case 3
                            % 125ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<125 
                                los3 = [los3; tr];
                            else 
                                win3 = [win3; tr];
                            end
                            
                        case 6
                            % 225ms is considered the threshold between distribution of reward times for win and lose trials
                             if eventTable.TotalRewardTime(tr)<225 
                                los6 = [los6; tr];
                             else 
                                win6 = [win6; tr];
                             end                           
                            
                        case 9
                            % 325ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<325 
                                los9 = [los9; tr];
                            else 
                                win9 = [win9; tr];
                            end
                    end
                end
            end
        end
        output = [output; struct('TrialIdx', los3, 'Value', 'lose 3', 'GroupingType', grType)];
        output = [output; struct('TrialIdx', los6, 'Value', 'lose 6', 'GroupingType', grType)];
        output = [output; struct('TrialIdx', los9, 'Value', 'lose 9', 'GroupingType', grType)];
        output = [output; struct('TrialIdx', win3, 'Value', 'win 3', 'GroupingType', grType)];
        output = [output; struct('TrialIdx', win6, 'Value', 'win 6', 'GroupingType', grType)];
        output = [output; struct('TrialIdx', win9, 'Value', 'win 9', 'GroupingType', grType)];
    
        
    case 'PreTrlOutcome & CurrentTrlVar'
        output = [];
        win0 = [];  win1 = [];   win4 = [];  
        los0 = [];  los1 = [];   los4 = [];  
        eventTable = struct2table(event);        
        for tr = 2:size(eventTable,1)  % for each trial
            if ~eventTable.TrialErrorCode(tr) % if current trial was successful
                if ( ~eventTable.TrialErrorCode(tr-1)  && ~strcmp(eventTable.type{tr-1}(1),'F') ) % if previous trial was successful (rewarded) && reawrd variance was not zero   
                    switch ['[' num2str(eventTable.RewardVariance(tr)) ' ' num2str(eventTable.expected_reward(tr)) ']']
                        case '[0 3]'
                            % 125ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<125 
                                los0 = [los0; tr];
                            else 
                                win0 = [win0; tr];
                            end
                        case '[0 6]'
                            % 225ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<225 
                                los0 = [los0; tr];
                            else 
                                win0 = [win0; tr];
                            end
                        case '[0 9]'
                            % 325ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<325 
                                los0 = [los0; tr];
                            else 
                                win0 = [win0; tr];
                            end
                        case '[1 3]'
                            % 125ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<125 
                                los1 = [los1; tr];
                            else 
                                win1 = [win1; tr];
                            end
                        case '[1 6]'
                            % 225ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<225 
                                los1 = [los1; tr];
                            else 
                                win1 = [win1; tr];
                            end
                        case '[1 9]'
                            % 325ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<325 
                                los1 = [los1; tr];
                            else 
                                win1 = [win1; tr];
                            end
                        case '[4 3]'
                            % 125ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<125 
                                los4 = [los4; tr];
                            else 
                                win4 = [win4; tr];
                            end
                        case '[4 6]'
                            % 225ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<225 
                                los4 = [los4; tr];
                            else 
                                win4 = [win4; tr];
                            end
                        case '[4 9]'
                            % 325ms is considered the threshold between distribution of reward times for win and lose trials
                            if eventTable.TotalRewardTime(tr)<325 
                                los4 = [los4; tr];
                            else 
                                win4 = [win4; tr];
                            end
                    end
                end
            end
        end
        output = [output; struct('TrialIdx', los0, 'Value', 'lose 0', 'GroupingType', grType)];
        output = [output; struct('TrialIdx', los1, 'Value', 'lose 1', 'GroupingType', grType)];
        output = [output; struct('TrialIdx', los4, 'Value', 'lose 4', 'GroupingType', grType)];
        output = [output; struct('TrialIdx', win0, 'Value', 'win 0', 'GroupingType', grType)];
        output = [output; struct('TrialIdx', win1, 'Value', 'win 1', 'GroupingType', grType)];
        output = [output; struct('TrialIdx', win4, 'Value', 'win 4', 'GroupingType', grType)];
    
    
    case 'PreTrialEV'
        ind1 = [];
        ind2 = [];
        ind3 = [];
        output = [];
        eventTable = struct2table(event);        
        for j = 2:size(eventTable,1)
            if ~eventTable.TrialErrorCode(j)
                if ~eventTable.TrialErrorCode(j-1)
                    if eventTable.expected_reward(j-1) == eventTable.expected_reward(j)
                        ind1 = [ind1;j];
                    elseif eventTable.expected_reward(j-1) > eventTable.expected_reward(j)
                        ind2 = [ind2;j];
                    elseif eventTable.expected_reward(j-1) < eventTable.expected_reward(j)
                        ind3 = [ind3;j];
                    end
                end
            end
        end

        output = [];
        output = [output; struct('TrialIdx',ind1,'Value', ['EV1 = EV2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind2,'Value', ['EV1 > EV2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind3,'Value', ['EV1 < EV2'],'GroupingType', grType)];
        
    case 'PreTrialVAR'
        ind1 = [];
        ind2 = [];
        ind3 = [];
        output = [];
        eventTable = struct2table(event);        
        for j = 2:size(eventTable,1)
            if ~eventTable.TrialErrorCode(j)
                if ~eventTable.TrialErrorCode(j-1)
                    if eventTable.RewardVariance(j-1) == eventTable.RewardVariance(j)
                        ind1 = [ind1;j];
                    elseif eventTable.RewardVariance(j-1) > eventTable.RewardVariance(j)
                        ind2 = [ind2;j];
                    elseif eventTable.RewardVariance(j-1) < eventTable.RewardVariance(j)
                        ind3 = [ind3;j];
                    end
                end
            end
        end

        output = [];
        output = [output; struct('TrialIdx',ind1,'Value', ['VAR1 = VAR2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind2,'Value', ['VAR1 > VAR2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind3,'Value', ['VAR1 < VAR2'],'GroupingType', grType)];
        
    case 'PreTrialEV&VAR'        
        ind = cell(1,9);
        output = [];
        eventTable = struct2table(event);        
        for j = 2:size(eventTable,1)
            if ~eventTable.TrialErrorCode(j)
                if ~eventTable.TrialErrorCode(j-1)
                    if eventTable.expected_reward(j-1) == eventTable.expected_reward(j)
                        if eventTable.RewardVariance(j-1) == eventTable.RewardVariance(j)
                            ind{1} = [ind{1};j];
                        elseif eventTable.RewardVariance(j-1) > eventTable.RewardVariance(j)
                            ind{2} = [ind{2};j];
                        elseif eventTable.RewardVariance(j-1) < eventTable.RewardVariance(j)
                            ind{3} = [ind{3};j];
                        end
                    elseif eventTable.expected_reward(j-1) > eventTable.expected_reward(j)
                        if eventTable.RewardVariance(j-1) == eventTable.RewardVariance(j)
                            ind{4} = [ind{4};j];
                        elseif eventTable.RewardVariance(j-1) > eventTable.RewardVariance(j)
                            ind{5} = [ind{5};j];
                        elseif eventTable.RewardVariance(j-1) < eventTable.RewardVariance(j)
                            ind{6} = [ind{6};j];
                        end
                    elseif eventTable.expected_reward(j-1) < eventTable.expected_reward(j)
                        if eventTable.RewardVariance(j-1) == eventTable.RewardVariance(j)
                            ind{7} = [ind{7};j];
                        elseif eventTable.RewardVariance(j-1) > eventTable.RewardVariance(j)
                            ind{8} = [ind{8};j];
                        elseif eventTable.RewardVariance(j-1) < eventTable.RewardVariance(j)
                            ind{9} = [ind{9};j];
                        end
                    end
                end
            end
        end

        output = [];
        output = [output; struct('TrialIdx',ind{1},'Value', ['EV1 = EV2  VAR1 = VAR2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind{2},'Value', ['EV1 = EV2  VAR1 > VAR2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind{3},'Value', ['EV1 = EV2  VAR1 < VAR2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind{4},'Value', ['EV1 > EV2  VAR1 = VAR2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind{5},'Value', ['EV1 > EV2  VAR1 > VAR2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind{6},'Value', ['EV1 > EV2  VAR1 < VAR2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind{7},'Value', ['EV1 < EV2  VAR1 = VAR2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind{8},'Value', ['EV1 < EV2  VAR1 > VAR2'],'GroupingType', grType)];
        output = [output; struct('TrialIdx',ind{9},'Value', ['EV1 < EV2  VAR1 < VAR2'],'GroupingType', grType)];
        
    case 'expecter_reward & RewardVariance'
        eventTable = struct2table(event);
        output = [];
        allConditions = unique([eventTable.expected_reward  eventTable.RewardVariance],'rows');
        % get the trial indices for variance 0
        idx = [];
        for gr = 1:size(allConditions,1)
            idx = find( (eventTable.expected_reward==allConditions(gr,1)) .* (eventTable.RewardVariance==allConditions(gr,2)) );
            output = [output; struct('TrialIdx', idx, ...
                                     'Value', allConditions(gr,:), ...
                                     'GroupingType', grType)];
            clear idx
       end
        clear gr
        
    case 'expecter_reward & RewardVariance & Outcome'
        eventTable = struct2table(event);
        output = [];
        allConditions = [3  0  0   %  -
                         6  0  0   % | there is no win or lose outcome for zero variance 
                         9  0  0   %  -
                         ...
                         3  1  -1  %  -
                         6  1  -1  % | variance 1; losing trials
                         9  1  -1  %  -
                         3  1  1   %  -
                         6  1  1   % | variance 1; winning trials
                         9  1  1   %  -
                         ...
                         3  4  -1  %  -
                         6  4  -1  % | variance 4; losing trials
                         9  4  -1  %  -
                         3  4  1   %  -
                         6  4  1   % | variance 4; winning trials
                         9  4  1]; %  -
        
         outcomeThreshold = [NaN; NaN; NaN; 125; 225; 325; 125; 225; 325; 125; 225; 325; 125; 225; 325];
                         
        % get the trial indices for variance 0
        rewdCueTrIdx = [];
        for gr = 1:size(allConditions,1)
            % find trials indeces for each group based on what the cue info
            rewdCueTrIdx = find( (eventTable.expected_reward==allConditions(gr,1)) .* (eventTable.RewardVariance==allConditions(gr,2)) );
            
            % treat conditions differently based on variablity of the reward 
            if ~isnan(outcomeThreshold(gr)) % this is for when reward variance is not zero  
                switch allConditions(gr,3)
                    case -1 % collect losing trials for this group
                        idx = rewdCueTrIdx( find(eventTable.TotalRewardTime(rewdCueTrIdx) <= outcomeThreshold(gr)) );
                    case 1  % collect winning trials for this group
                        idx = rewdCueTrIdx( find(eventTable.TotalRewardTime(rewdCueTrIdx) >= outcomeThreshold(gr)) );
                end
            else % this is for F3, F6 and F9 cue conditions
                idx = rewdCueTrIdx;
            end
            
            % pass this group to the output
            output = [output; struct('TrialIdx', idx, ...
                                     'Value', allConditions(gr,:), ...
                                     'GroupingType', grType)];

            clear rewdCueTrIdx idx
       end
       clear gr
       
    case 'PreTrlEV & PreTrlVar & CurrTrlEV' 
        % each row of the 'allConditions' contains following values -> [PreTrlEV   PreTrlVar  CurrTrlEV] 
        % in total there are 27 conditions
        allConditions = [3  0  3   %  -
                         3  0  6   % | 
                         3  0  9   %  -
                         3  1  3   %  -
                         3  1  6   % | 
                         3  1  9   %  -
                         3  4  3   %  -
                         3  4  6   % | 
                         3  4  9   %  -
                         ...
                         6  0  3   %  -
                         6  0  6   % | 
                         6  0  9   %  -
                         6  1  3   %  -
                         6  1  6   % | 
                         6  1  9   %  -
                         6  4  3   %  -
                         6  4  6   % | 
                         6  4  9   %  -                        
                         ...
                         9  0  3   %  -
                         9  0  6   % | 
                         9  0  9   %  -
                         9  1  3   %  -
                         9  1  6   % | 
                         9  1  9   %  -
                         9  4  3   %  -
                         9  4  6   % | 
                         9  4  9]; %  -                        
        
        output = [];
        eventTable = struct2table(event);
        idx = cell(27,1);
        for tr = 2:size(eventTable,1)
            if (~eventTable.TrialErrorCode(tr) && ~eventTable.TrialErrorCode(tr-1)) % if this trial and the previouse one were successful 
                currCondi = [eventTable.expected_reward(tr-1) eventTable.RewardVariance(tr-1) eventTable.expected_reward(tr)];    
                % find which condition/group it is
                cd = find( ismember(allConditions, currCondi, 'rows') );
                idx{cd} = [idx{cd}; tr];
                clear cd currCondi
            end
        end
        clear tr cd currCondi

        for gr = 1: length(idx)
            output = [output; struct('TrialIdx', idx{gr}, ...
                                     'Value', allConditions(gr,:), ...
                                     'GroupingType', grType)];
        end 
        
    case 'PreTrlEV & PreTrlOutcome & CurrTrlEV' 
        % each row of the 'allConditions' contains following values -> [PreTrlEV   PreTrlOutcome  CurrTrlEV] 
        % in total there are 27 conditions
        allConditions = [3  -1  3   %  -
                         3  -1  6   % | 
                         3  -1  9   %  -
                         3  0   3   %  -
                         3  0   6   % | 
                         3  0   9   %  -
                         3  1   3   %  -
                         3  1   6   % | 
                         3  1   9   %  -
                         ...
                         6  -1  3   %  -
                         6  -1  6   % | 
                         6  -1  9   %  -
                         6  0   3   %  -
                         6  0   6   % | 
                         6  0   9   %  -
                         6  1   3   %  -
                         6  1   6   % | 
                         6  1   9   %  -
                         ...
                         9  -1  3   %  -
                         9  -1  6   % | 
                         9  -1  9   %  -
                         9  0   3   %  -
                         9  0   6   % | 
                         9  0   9   %  -
                         9  1   3   %  -
                         9  1   6   % | 
                         9  1   9]; %  -                        
        
        output = [];
        eventTable = struct2table(event);
        idx = cell(27,1);
        for tr = 2:size(eventTable,1)
            if (~eventTable.TrialErrorCode(tr) && ~eventTable.TrialErrorCode(tr-1)) % if this trial and the previouse one were successful 
                currCondi = [eventTable.expected_reward(tr-1) GetTrialOutcome(eventTable(tr-1,:)) eventTable.expected_reward(tr)];    
                % find which condition/group it is
                cd = find( ismember(allConditions, currCondi, 'rows') );
                % collect current trial number for the relevant condition
                idx{cd} = [idx{cd}; tr];
                clear cd currCondi
            end
        end
        clear tr cd currCondi

        for gr = 1: length(idx)
            output = [output; struct('TrialIdx', idx{gr}, ...
                                     'Value', allConditions(gr,:), ...
                                     'GroupingType', grType)];
        end
        
        
        
    case 'PreTrlVar & PreTrlOutcome & CurrTrlEV' 
        % each row of the 'allConditions' contains following values -> [PreTrlEV   PreTrlOutcome  CurrTrlEV] 
        % in total there are 27 conditions
        allConditions = [0  0   3   %  -
                         0  0   6   % | 
                         0  0   9   %  -
                         ...
                         1  -1  3   %  -
                         1  -1  6   % | 
                         1  -1  9   %  -
                         1  1   3   %  -
                         1  1   6   % | 
                         1  1   9   %  -
                         ...
                         4  -1  3   %  -
                         4  -1  6   % | 
                         4  -1  9   %  -
                         4  1   3   %  -
                         4  1   6   % | 
                         4  1   9]; %  -                        
        
        output = [];
        eventTable = struct2table(event);
        idx = cell(size(allConditions, 1),1);
        for tr = 2:size(eventTable,1)
            if (~eventTable.TrialErrorCode(tr) && ~eventTable.TrialErrorCode(tr-1)) % if this trial and the previouse one were successful 
                currCondi = [eventTable.RewardVariance(tr-1) GetTrialOutcome(eventTable(tr-1,:)) eventTable.expected_reward(tr)];    
                % find which condition/group it is
                cd = find( ismember(allConditions, currCondi, 'rows') );
                % collect current trial number for the relevant condition
                idx{cd} = [idx{cd}; tr];
                clear cd currCondi
            end
        end
        clear tr cd currCondi

        for gr = 1: length(idx)
            output = [output; struct('TrialIdx', idx{gr}, ...
                                     'Value', allConditions(gr,:), ...
                                     'GroupingType', grType)];
        end 
        
        
    case 'CuePos & TargetPos & CurrEV & CurrVar' 
        % make a matrix with all posiible combinations of these variables
        allConditions = combvec([-1 1],[-1 1],[3 6 9],[0 1 4])';
        output = [];
        eventTable = struct2table(event);
        for gr = 1:size(allConditions,1)
            trlIdx = find((eventTable.cue_pos==allConditions(gr,1)) .* (eventTable.target_pos==allConditions(gr,2)) .* ...
                          (eventTable.expected_reward==allConditions(gr,3)) .* (eventTable.RewardVariance==allConditions(gr,4)));
            output = [output; struct('TrialIdx', trlIdx, ...
                                     'Value', allConditions(gr,:), ...
                                     'GroupingType', grType)];
            clear trlIdx                                 
        end 
        clear gr
        
    case 'CuePos & TargetPos & PreGood & CurrEV' 
        % make a matrix with all possible combinations of these variables
        allConditions = combvec([-1 1],[-1 1], [0,1] ,[0 0.25 0.5 0.75 1])';
        output = [];
        eventTable = struct2table(event);
        for gr = 1:size(allConditions,1)
            trlIdx = find((eventTable.cue_pos==allConditions(gr,1)) .* (eventTable.target_pos==allConditions(gr,2)) .* ...
                          (eventTable.pre_good==allConditions(gr,3)) .* (eventTable.expected_reward==allConditions(gr,4)));
            output = [output; struct('TrialIdx', trlIdx, ...
                                     'Value', allConditions(gr,:), ...
                                     'GroupingType', grType)];
            clear trlIdx                                 
        end 
        clear gr
      
    
    case 'CuePos & TargetPos & PreGood & CurrEV (uncertainty)' 
        % make a matrix with all possible combinations of these variables
        allConditions = combvec([-1 1],[-1 1], [0,1] ,[0 0.5 1])';
        output = [];
        eventTable = struct2table(event);
        for gr = 1:size(allConditions,1)
            trlIdx = find((eventTable.cue_pos==allConditions(gr,1)) .* (eventTable.target_pos==allConditions(gr,2)) .* ...
                          (eventTable.pre_good==allConditions(gr,3)) .* (eventTable.expected_reward==allConditions(gr,4)));
            output = [output; struct('TrialIdx', trlIdx, ...
                                     'Value', allConditions(gr,:), ...
                                     'GroupingType', grType)];
            clear trlIdx                                 
        end 
        clear gr
end

return