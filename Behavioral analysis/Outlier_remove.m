function [ totalTrials_new,outlier_indices,zscored_RT ] = Outlier_remove( totalTrials,outlierfield )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

zscore_thr = 2;
outlier_indices = [];

zscored_RT = zeros(length(totalTrials),1);

for i = 1:2  % '1' for Mojo and '2' for MacDuff
    monkey_indx = find([totalTrials.subjectID] == i);
    
    % remove outliers of event timings
    if isfield(outlierfield,'trial_start_time')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.startTimes;
        end
        outlier_indices1 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices1 = [];
    end
    clear data
    
    if isfield(outlierfield,'trial_fixation_on')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.fixTimes;
        end
        outlier_indices2 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices2 = [];
    end
    clear data
    
    if isfield(outlierfield,'trial_fixation_Acquired')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.fixAq;
        end
        outlier_indices3 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices3 = [];
    end
    clear data
    
    if isfield(outlierfield,'trial_cue_onset')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.cueTimes;
        end
        outlier_indices4 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices4 = [];
    end
    clear data
    
    if isfield(outlierfield,'trial_cue_offset')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.cueoffTimes;
        end
        outlier_indices5 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices5 = [];
    end
    clear data
    
    if isfield(outlierfield,'trial_target_onset')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.targTimes;
        end
        outlier_indices6 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices6 = [];
    end
    clear data    

    if isfield(outlierfield,'trial_target_acquired')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.targAqTimes;
        end
        outlier_indices7 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices7 = [];
    end
    clear data 
    
    if isfield(outlierfield,'trial_target_offset')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.targOffTimes;
        end
        outlier_indices8 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices8 = [];
    end
    clear data 
    
    if isfield(outlierfield,'trial_end_time')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.endTimes;
        end
        outlier_indices9 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices9 = [];
    end
    clear data 
    
    if isfield(outlierfield,'trial_reward_onset')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).RewardOnset;
        end
        outlier_indices10 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices10 = [];
    end
    clear data
    
    if isfield(outlierfield,'trial_reward_offset')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).RewardOffset;
        end
        outlier_indices11 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices11 = [];
    end
    clear data
    
    % remove outliers of inter-event timings
    if isfield(outlierfield,'trial_trial_length')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.endTimes - totalTrials(monkey_indx(j)).Time.startTimes;
        end
        outlier_indices12 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices12 = [];
    end
    clear data
    
    if isfield(outlierfield,'trial_fixation_to_acquired')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.fixAq - totalTrials(monkey_indx(j)).Time.fixTimes;
        end
        outlier_indices13 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices13 = [];
    end
    clear data    
  
    if isfield(outlierfield,'trial_cueonset_to_cueoffset')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).ActualEventTime.CueOffset - totalTrials(monkey_indx(j)).ActualEventTime.CueOnset;
        end
        outlier_indices14 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices14 = [];
    end
    clear data 
    
    if isfield(outlierfield,'trial_cueoffset_to_targetacquired')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.targAqTimes - totalTrials(monkey_indx(j)).Time.cueoffTimes;
        end
        outlier_indices15 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices15 = [];
    end
    clear data    
    
    if isfield(outlierfield,'trial_targeton_to_targetacquired')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).Time.targAqTimes - totalTrials(monkey_indx(j)).Time.targTimes;
        end
        outlier_indices16 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices16 = [];
    end
    clear data 
    
    if isfield(outlierfield,'trial_target_to_reward')
        for j = 1:numel(monkey_indx)
            data(j) = min(totalTrials(monkey_indx(j)).RewardOnTime) - totalTrials(monkey_indx(j)).ActualEventTime.TargetAcquired;
        end
        outlier_indices17 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices17 = [];
    end
    clear data
    
    if isfield(outlierfield,'trial_reward_duration')
        for j = 1:numel(monkey_indx)
            data(j) = totalTrials(monkey_indx(j)).RewardOffset - totalTrials(monkey_indx(j)).RewardOnset;
        end
        outlier_indices18 = monkey_indx(find(abs(zscore(data))>zscore_thr));
    else
        outlier_indices18 = [];
    end
    clear data
    % remove outliers of reaction time
    if isfield(outlierfield,'trial_reaction_time')
        [ z_scored ] = Z_score_trials( totalTrials(monkey_indx),'per_session','RT',0 );
        outlier_indices19 = monkey_indx(find(abs(z_scored)>zscore_thr));
        zscored_RT(monkey_indx) = z_scored;
    else
        outlier_indices19 = [];
    end

outlier_indices = [outlier_indices,unique([outlier_indices1,outlier_indices2,...
    outlier_indices3,outlier_indices4,outlier_indices5,outlier_indices6,...
    outlier_indices7,outlier_indices8,outlier_indices9,outlier_indices10,...
    outlier_indices11,outlier_indices12,outlier_indices13,outlier_indices14,...
    outlier_indices15,outlier_indices16,outlier_indices17,outlier_indices18,...
    outlier_indices19])];

clearvars -except outlier_indices zscore_thr totalTrials outlierfield zscored_RT
 
end

totalTrials(outlier_indices) = [];
totalTrials_new = totalTrials;

end

