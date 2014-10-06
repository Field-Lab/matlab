function [struct, success] = toStruct(data, e)
%TOSTRUCT toStruct method, returns data from electrode e
%   struct = toStruct(Data data, electrode) reads projection and spike
%   information 
%
%
%   tamachado@salk.edu 1/23/08


% Read in electrode e
try
    success = data.prjObject.readElectrode(e);
    
    if(success == false)
        data.badElectrodes(e) = true;
        error('Failure loading electrode %d', e);
    end

    % Save the data we just loaded into matlab
    data.prjStruct.spikeCount    =  data.prjObject.getSpikeCount();
    data.prjStruct.spikeTimes    =  data.prjObject.getSpikeTimes();
    data.prjStruct.projections   =  data.prjObject.getProjections();
    data.prjStruct.currElectrode =  data.prjObject.getElectrode();  
    
    % Kill trailing zeros in our data
    spikeCount = data.prjStruct.spikeCount;
    data.prjStruct.spikeTimes = data.prjStruct.spikeTimes(1:spikeCount);
    data.prjStruct.projections = data.prjStruct.projections(:,1:spikeCount);
     
    % Return the new data structure
    struct = data.prjStruct;

% Give up
catch
    error('Failure saving electrode %d', e);
end