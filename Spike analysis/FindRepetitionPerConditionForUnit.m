

% This function calls GroupTrials function, generates the output structure
% of the GroupTrials function; adds an additional field of binary values for each group which
% indicates for which trials of thar group, the unit was active.

function [output] = FindRepetitionPerConditionForUnit (groupingCriteria, event, unitPSTH)

[output] = GroupTrials(event, groupingCriteria);

% find for which trials of each group, the unit was active
for grp = 1:length(output)
    [output(grp).UnitActiveTrials , ~] = ismember( output(grp).TrialIdx , unitPSTH.active_eventStrctTrls );
end

return