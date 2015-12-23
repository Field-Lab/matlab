function [shortName, amp] = plotSelectivity(datarun, cellsToInclude, stimInfo, axesH)

% makes a mosaic plot with RF outlines colored according to response rate
%   inputs: 
%       datarun, after loading up STA fits (including whatever info is required by plot_rf_summaries.m)
%       cellsToInclude: a vector of cell IDs to be included in plot
%       stimInfo: structure containing info about electrical stimulation
%          stimInfo.pathToAnalysis - full path to where elecResp files are stored
%          stimInfo.suffix - elecResp suffix, typically indicating pulse phase width
%          stimInfo.stimElec -  electrode number of stimulus (typically this
%                               is the same as the pattern number for
%                               1-elec stim)
%          stimInfo.movieNo - movie number of stimulus

plotColors = jet(101);

pathToAnalysis = stimInfo.pathToAnalysis;
stimElec = stimInfo.stimElec;
movie = stimInfo.movieNo;
cell = stimInfo.id;
if isfield(stimInfo, 'suffix')
    suffix = stimInfo.suffix;
else
    suffix = [];
end

if ~isempty(axesH)
    axes(axesH)
else
    figure
end

%figH = figure
hold on
for ii = 1:length(cellsToInclude)
    load([pathToAnalysis 'elecResp_n' num2str(cellsToInclude(ii)) '_p' num2str(stimElec) suffix])

    movieInd = find(elecResp.stimInfo.movieNos == movie);
    succRate = round(100*sum(elecResp.analysis.latencies{movieInd}~=0)/length(elecResp.analysis.latencies{movieInd}));
    if cellsToInclude(ii) == cell;
        shortName = elecResp.names.rrs_short_name;
        amp = elecResp.stimInfo.stimAmps(movieInd);
    end
    clear elecResp
    
%     if succRate == 100 % 100% success rate is included in final color
%         plot_rf_summaries(datarun, cellsToInclude(ii), 'fit_width', 3, 'fit_color', plotColors(succRate+1,:))
%     else
    plot_rf_summaries(datarun, cellsToInclude(ii), 'fit_width', 3, 'fit_color', plotColors(succRate+1,:))
%     end
end
hold off





