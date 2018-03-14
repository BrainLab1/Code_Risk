
% This function returns spike counts in the specified time interval around the specified event

% grpIn:        has the format of output of the GroupTrials() output
% spkTimeData:  spike times relative to the alignedEvent
% timeInterval: time interval in seconds relative to the alignedEvent to count the number of spikes  


function [output] = GetSpikeCount (grpIn, spkTimeData, alignedEvent, timeInterval)

% generate the spike count field name based on the event name and time
% % % % spkTimeFieldName = [ 'spkCount_' alignedEvent '_' num2str(timeInterval(1)) '-' num2str(timeInterval(2)) ];


for grp = 1:length(grpIn)
    aa = ismember( spkTimeData.active_eventStrctTrls , grpIn(grp).TrialIdx(grpIn(grp).UnitActiveTrials) ); 
    
    % initialize spike times with NaN values
    grpIn(grp).AlignedSpikeTimes = mat2cell( nan(length(grpIn(grp).TrialIdx),1), ones(length(grpIn(grp).TrialIdx),1), 1 );
    
    grpIn(grp).AlignedSpikeTimes( grpIn(grp).UnitActiveTrials ) = spkTimeData.SpikeTimeByTrial(aa);
    
    % initialize spike counts with NaN values
    grpIn(grp).SpikeCount = nan(length(grpIn(grp).TrialIdx),1);
    
    % count spikes in the interval
    for tr = 1:length(grpIn(grp).AlignedSpikeTimes)
        
        if isnan(grpIn(grp).AlignedSpikeTimes{tr})
            grpIn(grp).SpikeCount(tr,1) = NaN;
        else if isempty(grpIn(grp).AlignedSpikeTimes{tr})
                grpIn(grp).SpikeCount(tr,1) = 0;
            else
                a1 = grpIn(grp).AlignedSpikeTimes{tr}>= timeInterval(1);
                a2 = grpIn(grp).AlignedSpikeTimes{tr}<= timeInterval(2);
                grpIn(grp).SpikeCount(tr,1) =  sum(a1 .* a2);
                clear a1 a2
            end
        end
        
    end
    clear tr aa
    
end
clear grp

output = grpIn;
return

