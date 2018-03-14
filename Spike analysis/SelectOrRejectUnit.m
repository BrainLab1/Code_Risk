

% crreated on 12.02.2018

% select or reject the unit based on the criteria indicated in the function input
% Criteria:
%    1 - 'RepetitionPerCondition'
%        - dataIn: should be the output of the  FindRepetitionPerConditionForUnit() function  
%        - options: is a structure with field name 'MinRepPerCnd' which shows minimum repetition per condition the unit must have to be selected 
%        - fncOut: integer variable which takes values 1 (if neuron passed the criteria) or 0 (if neuron did not pass the criteria)



function [fncOut, decision] = SelectOrRejectUnit(dataIn, criteria, options)

switch criteria
    
    case 'RepetitionPerCondition'
        repPerCnd = nan(length(dataIn),1);
        repPerCnd = cellfun(@(x) sum(x) , {dataIn.UnitActiveTrials}');
        fncOut = dataIn;
        
        % add number of repetitions per condition/group to the input data and pass it to the output
        for grp = 1:length(dataIn)
            fncOut(grp).RepPerCnd = repPerCnd(grp);
        end
        clear grp
        
        % report whether or not the unit passed the criteria 
        if sum(repPerCnd >= options.MinRepPerCnd) == length(dataIn)
           decision = 1;
        else 
           decision = 0;
        end 
           
        clear repPerCnd
end

return