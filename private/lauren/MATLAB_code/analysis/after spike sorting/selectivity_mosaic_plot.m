function selectivity_mosaic_plot(datarun, goodCells, stimElec, movieNo, pathToAnalysis, varargin)

% datarun: standard datarun structure containing the fields necessary for plot_rf_summaries:
%      .ei.position
%      .piece.T_monitor_to_array
%
% goodCells: cell array of vectors of cell IDs in each cell class (ON par, OFF par, ON
% midg, OFF midg) to be plotted, corresponding with datarun
%
% stimElec: electrode used for stimulation
% movieNo: movie number corresponding to chosen stimulus amplitude
% pathToAnalysis: full path to the directory containing the elecResp files
% elecRespSuffix: if there is anything to append to elecResp filename after _p## (e.g. _w50)
%
%

p = inputParser;

p.addRequired('datarun', @isstruct)
p.addRequired('goodCells', @iscell)
p.addRequired('stimElec', @isnumeric)
p.addRequired('movieNo', @isnumeric)
p.addRequired('pathToAnalysis', @ischar)

p.addParamValue('elecRespSuffix', [], @ischar)
p.addParamValue('succRateOverride', [], @isnumeric) % allows user to provide success rate value rather than get it from elecResp file
% must be in the format: [cellID1 succRate1; cellID2 succRate2;...]
p.addParamValue('plotElecs', false, @islogical)

p.parse(datarun, goodCells, stimElec, movieNo, pathToAnalysis, varargin{:})

elecRespSuffix = p.Results.elecRespSuffix;
succRateOR = p.Results.succRateOverride;
plotElecs = p.Results.plotElecs;




%determine electrode positions in STA space
T_array_to_monitor = fliptform(datarun.piece.T_monitor_to_array);
T_monitor_to_STA = fliptform(coordinate_transform(datarun, 'monitor'));
T = maketform('composite',[T_monitor_to_STA T_array_to_monitor]);
e_pos = tformfwd(T, datarun.ei.position);

%
axesPositions = {[0 0 0.5 0.5], [0.5 0 0.5 0.5], [0 0.5 0.5 0.5], [0.5 0.5 0.5 0.5]};
plotColors = jet(101);


%mosaics plot
figure('position', [100 100 300 300])
for ii = 1:4
    axes('position', axesPositions{ii})
    hold on
    for jj = 1:length(goodCells{ii})
        if isempty(succRateOR) || ~any(succRateOR(:,1)==goodCells{ii}(jj))

            load([pathToAnalysis 'elecResp_n' num2str(goodCells{ii}(jj)) '_p' num2str(stimElec) elecRespSuffix])
            movieInd = find(elecResp.stimInfo.movieNos == movieNo);

            if ~elecResp.analysis.finalized(movieInd)
                warning(['movie ' num2str(movieNo) ' not finalized for elecResp_n' num2str(goodCells{ii}(jj)) '_p' num2str(stimElec) elecRespSuffix])
            end
            succRate = round(100*sum(elecResp.analysis.latencies{movieInd}~=0)/length(elecResp.analysis.latencies{movieInd}));

            clear elecResp
        else %user-specified success rate
            succRate = round(100*succRateOR(succRateOR(:,1)==goodCells{ii}(jj),2));
        end
        if jj == 1 && ii == 1 %include array outline
            plot_rf_summaries(datarun, goodCells{ii}(jj), 'fit_width', 1, 'fit_color', [0 0 0], 'fill_color', plotColors(succRate+1,:), 'array', true)
        else
            plot_rf_summaries(datarun, goodCells{ii}(jj), 'fit_width', 1, 'fit_color', [0 0 0], 'fill_color', plotColors(succRate+1,:), 'array', true)
        end
    end
    
    if plotElecs
        plot(e_pos(:,1), e_pos(:,2), 'k.')
        plot(e_pos(stimElec,1), e_pos(stimElec,2), 'r*')
    end
    
    hold off
    axis equal
    set(gca, 'xlim', [0 30], 'ylim', [0 30])
    axis off
end


