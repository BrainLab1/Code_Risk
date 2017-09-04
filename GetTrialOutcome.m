
% created on 03.09.2017 by Bahareh

% This function returns the trial outcome (win or loose) for one trial based on the
% total reward time and cue type. the input has to be a structure eith fileds
% .type (i.e., cue types of the risk task) and .TotalRewardTime
% output will be 1, 0 or -1 indicating win, neutal (for fixed reward trials) or loose outcome, respectively

function [outcome] = GetTrialOutcome(in)
outcome = NaN;
meanThreshold = [125, 225, 325]; % these are approximate reward time for fixed reward trials of three expected values   

% get the expected reawrd from the cue type and the first character of cue type identifier 
if iscell(in.type)
    expRwd = str2num(in.type{1}(2));
    cueType = in.type{1}(1);
else
    expRwd = str2num(in.type(2));
    cueType = in.type(1);
end

idx = expRwd/3; % idx will be 1, 2 or 3 indicating expected values of 3, 6 or 9

if (strcmp(cueType,'F')) % zero-variance trials
    outcome = 0;
else if in.TotalRewardTime < meanThreshold(idx) % loosing trials
        outcome = -1;
    else  % win trials
        outcome = 1;
    end
end

return