function [analyzedCells,activationProb] = plotActivatedCells(fPath,eipath,patternNo,varargin)
% PLOTACTIVATEDCELLS() plots a group of cells that are activated at a given
% stimulation electrode, at the stimulation amplitude that activates one
% cell with 0.75 probability
% inputs:       fPath
%               eipath
%               patternNo
% optional:     probability number between 0 and 1
%               cellIDs list to analyze just a subset
%               plotsOff  turn plots off; 

% Set up defaults for optional parameters
probability = 0.75;
cellIDs = 1:1000; % Dummy list so that all possible values will be included
plotsOff = 0; 

nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

% Read the optional input arguments
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'probability'
            probability = varargin{j*2};
        case 'cellids'
            cellIDs = varargin{j*2};
        case 'plotsoff'
            plotsOff = varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

dirInfo = dir(fPath);
if ~plotsOff
    fh =figure;
    fh.Position =  [114         612        1534         486];
    ax1 = subplot(1,2,1);
    ax2 = subplot(1,2,2);
    ax1.XLabel.String = 'stimulation amplitude';
    ax1.YLabel.String = 'activation probability';
    colors = lines(length(cellIDs));
end
[eiM,nIdList] = convertEiFileToMatrix(eipath);
% could potentially change to just specify the elecRespAuto filename to
% load, and wouldn't have to search through an entire (large) directory to
% find the correct file names. But this is only good for when we specify
% the cellIDs. 
nn = 1; % increment the neuron that has been plotted. 
ph_all=[];
celllist = {};
allFits = []; 
for  n = 3:size(dirInfo,1) 
    if ~dirInfo(n).isdir 
        fname = dirInfo(n).name;
        i = find(fname=='_',2,'first');       
        if strcmp(['p' num2str(patternNo)],fname(i(2)+1:end-4)) && ismember(str2double(fname(i(1)+2:i(2)-1)),cellIDs)
            temp = load([fPath fname]);
            elecResp = temp.elecRespAuto; clear temp;
            try
                spikes = cell2mat(elecResp.spikes);
                responseProb = sum(spikes,2)/size(spikes,2);
            catch err
                disp([fPath fname ' is empty']);
                throw(err)
            end

            [~, completeFit, ~] = fitToErf_inline(elecResp,...
                responseProb,0);
            if ~plotsOff
                axes(ax1); hold on; ph = plot(completeFit(1,:),completeFit(2,:),'LineWidth',2,'Color',colors(nn,:));
                ph_all = [ph_all ph]; %#ok<*AGROW>
            end
            allFits{nn} = completeFit; 
            celllist = [celllist; num2str(elecResp.neuronInfo.neuronIds)];

            idx = find(completeFit(2,:)>probability,1,'first');
            if isempty(idx)
                stimAmpForGivenProb(nn) = 5; 
            else
                stimAmpForGivenProb(nn) = completeFit(1,idx); 
            end
            nn = nn+1;
        end % end load elecRespAuto files from a particular neuron
    end
end % End search through directory
if ~plotsOff
    axes(ax1); line(ax1.XLim,probability*[1 1],'LineStyle','--','Color','black');
    legend(ph_all,celllist);
    axes(ax2);
end
try
    stimAmpToPlot = min(stimAmpForGivenProb);
catch
    stimAmpToPlot = 0;
end
if ~plotsOff
    colors = jet(10);
    for ii = 1:size(allFits,2)
        currentFit = allFits{ii};
        currentProb(ii) = currentFit(2,find(currentFit(1,:) == stimAmpToPlot));
        colorIdx = ceil(currentProb(ii)*10);
        if colorIdx == 0; colorIdx = 1; end
        eiAmps = max(eiM(:,:,find(nIdList == str2double(celllist(ii)))),[],2);
        [axon_x, axon_y, ~, soma_x, soma_y] = weighted_axon_poly_reg(eiAmps);
        plot(axon_x, axon_y, 'Color', colors(colorIdx,:),'LineWidth',2, 'DisplayName', celllist{ii});
        legend('-DynamicLegend'); hold on;
        h = scatter(soma_x,soma_y,max(eiAmps)+0.1,colors(colorIdx,:),'filled');
        set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    colorbar;
    axes(ax2);
    axis off; axis image;
    [xc,yc] = getElectrodeCoords512();
    stimElec = unique(elecResp.stimInfo.listStimElecs);
    scatter(xc(stimElec),yc(stimElec),50,'k','filled');
    line(stimAmpToPlot*[1 1],ax1.YLim,'LineStyle','--','Color','black');
else
    try
    for ii = 1:size(allFits,2)
        
        try
            currentFit = allFits{ii};
        currentProb(ii) = currentFit(2,find(currentFit(1,:) == stimAmpToPlot));
        catch
            currentProb(ii) =0;
        end
    end
    catch
        currentProb = 0; 
        disp(['warning - no cells for this pattern' num2str(patternNo)])
    end
end
analyzedCells = str2double(celllist); 
try
    activationProb = currentProb;
catch
    activationProb = 0;
end
    function [threshold,completeFit, erfErr] = fitToErf_inline(elecResp,responseProb,plotResponseCurves)
        nMovies = length(elecResp.stimInfo.listAmps);
        stimIdx = find(abs(elecResp.stimInfo.listAmps(1,:)) == max(abs(elecResp.stimInfo.listAmps(1,:))));
        data = zeros(2, nMovies);
        data(1,:) = elecResp.stimInfo.listAmps(:,stimIdx); %#ok<FNDSB>
        data(2,:) = responseProb;
        data(3,:) = elecResp.tracesInfo.I;
               
        % linear-based
        data(1,:) = abs(data(1,:));
        try
            [erfParams, completeFit, erfErr] = erfFitter(data, 2, -1, 'makePlot', plotResponseCurves);
            threshold = -erfParams(2)/erfParams(1);
        catch
            threshold = 0;
            completeFit = zeros(2,size(data,2)); 
            erfErr = 0;
        end
        
    end
end% end function