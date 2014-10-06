function dsInfo = badElecCheck(varargin)

%function to check for bad recording electrodes

p = inputParser;



p.addParamValue('badElecThresh', 0.5, @isnumeric)
p.addParamValue('makePlots', true, @islogical)

p.parse(varargin{:})

params = p.Results;


% badElecThresh: fraction of mean(over electrodes) summed EI amplitude that an electrode
% must fall below to be considered "bad"

dsInfo = dataset_list_axon_directions_StanfordDirs();

[xCoords yCoords] = getElectrodeCoords61();

for ii = 1:length(dsInfo)
    eiInd = strfind(dsInfo(ii).eiPath, '.ei');
    datarun = load_data(dsInfo(ii).eiPath(1:eiInd-1));
    datarun = load_neurons(datarun);
    datarun = load_ei(datarun, 'all');
    
    eiAmpSum = zeros(64,1);
    for jj = 1:length(datarun.ei.eis)
        eiAmpSum = eiAmpSum + max(abs(datarun.ei.eis{jj}),[],2);
    end
    
    disconBin = false(64,1);
    disconBin([9 25 57]) = true;
    dsInfo(ii).badElecs = find(eiAmpSum < params.badElecThresh*mean(eiAmpSum(~disconBin)) & ~disconBin);
    
    if params.makePlots
        figure; hold on
        title(dsInfo(ii).eiPath)
        for jj = 1:64
            if ~any([9 25 57]==jj)
                if ~any(dsInfo(ii).badElecs == jj)
                    plot(xCoords(jj), yCoords(jj), 'o', 'markersize', 50*eiAmpSum(jj)/max(eiAmpSum), 'markerEdgeColor', [0.5 0.5 0.5])
                else
                    plot(xCoords(jj), yCoords(jj), 'ro', 'markersize', 50*eiAmpSum(jj)/max(eiAmpSum))
                end
                text(xCoords(jj), yCoords(jj), num2str(jj))
            end
        end
        axis equal
    end
    %keyboard
end