function [ z_scored ] = Z_score_trials( totalTrials,type,field_name,option )
%UNTITLED2 Summary of this function goes here
%   session based z-score calculation
%   option.alignment
%   option.interval
%   totalTrials contains all completed trials

% type: 'per_session', 'whole'


if strcmp(type,'per_session')

[output] = GroupTrials(totalTrials, 'SessionID');

switch field_name
    case 'RT'
        z_scored = zeros(length(totalTrials),1);
        for i = 1:numel(output)
            ind = output(i).TrialIdx;
            tmp = zscore(real(log([totalTrials(ind).DiodeReactionTime])));
            z_scored(ind) = tmp;
%             z_scored = [z_scored,tmp];
%             tmp2(ind) = tmp;
        end
%         z_scored = tmp2;
    case 'lickData'
        switch option.alignment
            case 'Cue'
                % cue aligned
                for i = 1:length(totalTrials)
                    cue_ind = totalTrials(i).Time.cueTimes;
                    lick_sig = totalTrials(i).lickData;
                    aligned_lick(i,:) = lick_sig(cue_ind + 1000*min(option.interval) : cue_ind + 1000*max(option.interval));
                    clear lick_sig
                end
            case 'Reward'
                % reward aligned
                for i = 1:length(totalTrials)
                    reward_ind = totalTrials(i).RewardOnset;
                    lick_sig = totalTrials(i).lickData;
                    aligned_lick(i,:) = lick_sig(reward_ind + 1000*min(option.interval) : reward_ind + 1000*max(option.interval));
                    clear lick_sig
                end
        end
        z_scored = [];
        for i = 1:numel(output)
            ind = output(i).TrialIdx;
            z_scored = [z_scored;zscore(aligned_lick(ind,:))];
        end
end


else
 
switch field_name
    case 'RT'
        z_scored = zscore([totalTrials.DiodeReactionTime]);
    case 'lickData'
        switch option.alignment
            case 'Cue'
                % cue aligned
                for i = 1:length(totalTrials)
                    cue_ind = totalTrials(i).Time.cueTimes;
                    lick_sig = totalTrials(i).lickData;
                    aligned_lick(i,:) = lick_sig(cue_ind + 1000*min(option.interval) : cue_ind + 1000*max(option.interval));
                    clear lick_sig
                end
            case 'Reward'
                % reward aligned
                for i = 1:length(totalTrials)
                    reward_ind = totalTrials(i).RewardOnset;
                    lick_sig = totalTrials(i).lickData;
                    aligned_lick(i,:) = lick_sig(reward_ind + 1000*min(option.interval) : reward_ind + 1000*max(option.interval));
                    clear lick_sig
                end
        end
        z_scored = zscore(aligned_lick);
end   
    
end