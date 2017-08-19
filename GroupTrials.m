
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
        
    case {'value', 'expected_reward', 'cue_pos', 'target_pos', 'pre_good', 'taskCondCode','subjectID','SessionID'}
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
        
        
end

return