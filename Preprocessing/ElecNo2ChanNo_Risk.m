

% last update: 05.02.2018 by Bahareh
% this function converts electrode number to channel number for the risk
% task; it will be used in CutAndAlignSpikes

function chNum = ElecNo2ChanNo_Risk(Monkey, elecNo)

        % get the electrode array map
        eval(['[elecArray1, elecArray2] = electrodepinout' Monkey '(''elec'');']) 
        % get the channel array map
        eval(['[chArray1, chArray2] = electrodepinout' Monkey '(''chan'');']) 
        
        % check if the electrode belongs to Array1
        [row,col] = find(elecNo==elecArray1);
        if ~isempty(row)  % if the electrode belongs to array 1
            chNum = chArray1(row,col);
        else % if the electrode belongs to array 2
            [row,col] = find(elecNo==elecArray2);
            chNum = chArray2(row,col);     
        end

      
return