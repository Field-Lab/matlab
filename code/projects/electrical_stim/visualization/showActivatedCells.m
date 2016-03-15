function [activatedCells,activationThresh,nonActivatedCells] = showActivatedCells(fPath,electrodeNo)

dirInfo = dir(fPath);
activatedCells = [];
activationThresh = []; 
nonActivatedCells = []; 

for  n = 3:size(dirInfo,1) 
    if ~dirInfo(n).isdir 
        fname = dirInfo(n).name;
        i = find(fname=='_',2,'first');
        if strcmp(['p' num2str(electrodeNo)],fname(i(2)+1:end-4))
            temp = load([fPath fname]);
            elecResp = temp.elecRespAuto; clear temp;
            try
                spikes = cell2mat(elecResp.spikes);
                responseProb = sum(spikes,2)/size(spikes,2);
%                 stimAmps     = abs(elecResp.stimInfo.listAmps);
            catch err
                disp([fPath fname ' is empty']);
                throw(err)
            end

            [threshold, completeFit, erfErr] = fitToErf_inline(elecResp,...
                responseProb,0);
            if threshold > 0 && threshold < 5
                activatedCells = [activatedCells elecResp.neuronInfo.neuronIds];
                activationThresh = [activationThresh threshold]; 
            else % Case where pattern did not activate the cell
                nonActivatedCells = [nonActivatedCells elecResp.neuronInfo.neuronIds];
            end
            
        end % end load elecRespAuto files from a particular neuron
    end
end % End search through directory
%
    function [threshold,completeFit, erfErr] = fitToErf_inline(elecResp,responseProb,plotResponseCurves)
        nMovies = length(elecResp.stimInfo.listAmps);
        
        data = zeros(2, nMovies);
        data(1,:) = elecResp.stimInfo.listAmps;
        data(2,:) = responseProb;
        data(3,:) = elecResp.tracesInfo.I;
        
        
        % linear-based
        data(1,:) = abs(data(1,:));
        try
            [erfParams, completeFit, erfErr] = erfFitter(data, 2, -1, 'makePlot', plotResponseCurves);
            threshold = -erfParams(2)/erfParams(1);
        catch
            threshold = 0;
            completeFit = zeros(1,size(data,2)); 
            erfErr = 0;
        end
        
    end
end% end function