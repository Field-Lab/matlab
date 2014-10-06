function cell_types_overview_fig()
%%only a function so that it can have nested functions

clear all

figure('position', [100 100 1060 750])

a.raw(1) = axes('units', 'pixels', 'position', [150 585 60 45]);
a.raw(2) = axes('units', 'pixels', 'position', [350 585 60 45]);
a.raw(3) = axes('units', 'pixels', 'position', [550 585 60 45]);
a.raw(4) = axes('units', 'pixels', 'position', [750 585 60 45]);
a.raw(5) = axes('units', 'pixels', 'position', [950 585 60 45]);

a.sub(1) = axes('units', 'pixels', 'position', [60  580 160 120]);
a.sub(2) = axes('units', 'pixels', 'position', [260 580 160 120]);
a.sub(3) = axes('units', 'pixels', 'position', [460 580 160 120]);
a.sub(4) = axes('units', 'pixels', 'position', [660 580 160 120]);
a.sub(5) = axes('units', 'pixels', 'position', [860 580 160 120]);

a.raster(1) = axes('units', 'pixels', 'position', [60  420 160 120]);
a.raster(2) = axes('units', 'pixels', 'position', [260 420 160 120]);
a.raster(3) = axes('units', 'pixels', 'position', [460 420 160 120]);
a.raster(4) = axes('units', 'pixels', 'position', [660 420 160 120]);
a.raster(5) = axes('units', 'pixels', 'position', [860 420 160 120]);

a.curve(1) = axes('units', 'pixels', 'position', [60  240 160 140]);
a.curve(2) = axes('units', 'pixels', 'position', [260 240 160 140]);
a.curve(3) = axes('units', 'pixels', 'position', [460 240 160 140]);
a.curve(4) = axes('units', 'pixels', 'position', [660 240 160 140]);
a.curve(5) = axes('units', 'pixels', 'position', [860 240 160 140]);

a.sel(1) = axes('units', 'pixels', 'position', [60  30 160 160]);
a.sel(2) = axes('units', 'pixels', 'position', [260 30 160 160]);
a.sel(3) = axes('units', 'pixels', 'position', [460 30 160 160]);
a.sel(4) = axes('units', 'pixels', 'position', [660 30 160 160]);
a.sel(5) = axes('units', 'pixels', 'position', [860 30 160 160]);

a.colorBar   = axes('units', 'pixels', 'position', [10 30 40 160]);

rasterOptions.markerSize = 5; %for raster plot
rasterOptions.markerColor = [0 0 0];

mosaicScaleXLim{1} = [6 26];
mosaicScaleYLim{1} = [9 29];
mosaicScaleXLim{2} = [5.5 25.5];
mosaicScaleYLim{2} = [5.5 25.5];
mosaicScaleXLim{3} = [2.5 22.5];
mosaicScaleYLim{3} = [6.5 26.5];
mosaicScaleXLim{4} = [4 24];
mosaicScaleYLim{4} = [10 30];
mosaicScaleXLim{5} = [5 25];
mosaicScaleYLim{5} = [3 23];

selMovNos = [50 52 22 39 47]; %to be used in selectivity plots
tracesMovNos = [44 48 17 31 43]; %to be used in traces and raster plots

respCurveXLim = {[0 150], [0 150], [0 150], [0 100], [0 150]};
%eiRegions = {7:46, 9:48, 7:46, 13:52, 9:48};
eiRegions = {8:48, 10:50, 8:48, 14:54, 10:50};
subYLim = {[-350 200], [-175 100], [-87.5 50], [-175 100], [-350 200]};
subYTick = {[-200 0 200], [-100 0 100], [-50 0 50], [-100 0 100], [-200 0 200]};

rawYLim = {[-1140 230], [-1080 -460], [-890 -120], [-1020 -450], [-1320 -340]};

cellNames = {'ON parasol', 'OFF parasol', 'ON midget', 'OFF midget', 'SBC'};

[cellInfoAll datarunAll] = cell_list_overview_fig();

%%
for ii = 1:5
    if ii ~= 5
        datarun = load_data(datarunAll(ii).pathToData);
        datarun = load_index(datarun);
        datarun = load_neurons(datarun);
        datarun = load_params(datarun, struct('verbose',1));
        datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
        
        datarun.names.obvius_fit_path = datarunAll(ii).obvius_fit_path;
        datarun.cellsToInclude = datarunAll(ii).cellsToInclude;
        
        %use obvius fits because they're nicer!
        datarun = load_obvius_sta_fits(datarun);
        datarun = get_sta_fits_from_obvius(datarun, datarun.cellsToInclude);
        
    else %SBC: special case
        %load up blurred and unblurred fits in order to figure out scaling of RFs
        
        %fits to unblurred STAs
        datarun_unblurred = load_data(datarunAll(5).pathToData);
        
        datarun_unblurred = load_index(datarun_unblurred);
        datarun_unblurred = load_neurons(datarun_unblurred);
        datarun_unblurred = load_params(datarun_unblurred, struct('verbose',1));
        datarun_unblurred = load_sta(datarun_unblurred,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
        
        datarun_unblurred.names.obvius_fit_path = datarunAll(5).obvius_fit_path;
        datarun_unblurred.cellsToInclude = datarunAll(5).cellsToInclude;
        
        datarun_unblurred = load_obvius_sta_fits(datarun_unblurred);
        datarun_unblurred = get_sta_fits_from_obvius(datarun_unblurred, datarun_unblurred.cellsToInclude);
        
        original_fits = [datarun_unblurred.obvius.sta_fits{get_cell_indices(datarun_unblurred, datarunAll(5).matchedFits)}];
        original_fits = [original_fits.sd]; %returns vector of SD pairs: [cell1SD1 cell1SD2 cell2SD1 cell2SD2...]
        original_mean_area = pi*mean(original_fits(1:2:end).*original_fits(2:2:end));% area of ellipes = pi*SD1*SD2
        
        %fits to blurred STAs
        datarun = load_data(datarunAll(5).pathToDataBlur);
        datarun = load_index(datarun);
        datarun = load_neurons(datarun);
        datarun = load_params(datarun, struct('verbose',1));
        datarun = load_sta(datarun,'load_sta',[],'save_sta',0,'save_rf',1,'verbose',1);
                
        datarun.names.obvius_fit_path = datarunAll(5).obvius_fit_path_blur;
        datarun.cellsToInclude = datarunAll(5).cellsToInclude;
        
        datarun = load_obvius_sta_fits(datarun);
        datarun = get_sta_fits_from_obvius(datarun, datarun.cellsToInclude);
        
        blurred_fits = [datarun.obvius.sta_fits{get_cell_indices(datarun, datarunAll(5).matchedFits)}];
        blurred_fits = [blurred_fits.sd];
        blurred_mean_area = pi*mean(blurred_fits(1:2:end).*blurred_fits(2:2:end));
        
        scale_factor = sqrt(original_mean_area/blurred_mean_area);
        
        %tweak blurred fit diameters to have equal area to unblurred version (based
        %on cells listed in matchedFits)
        for jj = 1:length(datarun.stas.fits)
            if ~isempty(datarun.stas.fits{jj})
                datarun.obvius.sta_fits{jj}.sd = scale_factor*datarun.obvius.sta_fits{jj}.sd;
                datarun.stas.fits{jj}.sd = scale_factor*datarun.stas.fits{jj}.sd;
            end
        end
    end
        
    cellInfo = cellInfoAll(ii);
    cellInfo.movieNo = selMovNos(ii);
    tracesMovie = tracesMovNos(ii); %to be used in traces and raster plots
    
    %selectivity plot
    plotSelectivity(datarun, datarun.cellsToInclude, cellInfo, a.sel(ii));
    set(a.sel(ii), 'xTick', [], 'yTick', [], 'box', 'on', 'xlim', mosaicScaleXLim{ii}, 'ylim', mosaicScaleYLim{ii})
        
    %load elecResp for rest of the plots
    elecResp = []; %because it's a static workspace
    load([cellInfo.pathToAnalysis 'elecResp_n' num2str(cellInfo.id) '_p' num2str(cellInfo.stimElec) cellInfo.suffix]);
    
    %response curve plot
    plotRespCurve(elecResp, a.curve(ii), [cellInfo.movieNo tracesMovie], cellInfo.PW)
    set(a.curve(ii), 'xlim', respCurveXLim{ii}, 'ylim', [0 1])
    if ii==1; ylabel('response probability'); end
    
    %raster plot
    if ii == 1
        rasterOptions.labelYAxis = true;
    else
        rasterOptions.labelYAxis = false;
    end
    
    %subtract 2 from nonzero latency values because stimulus is applied at sample 2
    elecRespCopy = elecResp; %so that it isn't accidentally saved out with altered latency values
    for jj = 1:length(elecResp.stimInfo.movieNos)
        if ~isempty(elecResp.analysis.latencies{jj})
            elecRespCopy.analysis.latencies{jj}(elecResp.analysis.latencies{jj}~=0) =...
                elecResp.analysis.latencies{jj}(elecResp.analysis.latencies{jj}~=0)-2;
        end
    end
    
    plot_raster(a.raster(ii), elecRespCopy, tracesMovie, rasterOptions)
    set(gca, 'xLim', [0 2], 'xtick', [0 1 2])
        
    %plots of traces
    plotTraces(elecResp, false,  tracesMovie, cellInfo.gain, a.sub(ii))
    hold on; plot(0:0.05:2, (1000/cellInfo.gain)*elecResp.cells.mainEI(elecResp.cells.recElec,eiRegions{ii}), 'k--', 'linewidth', 2); hold off
    
    title(a.sub(ii), cellNames{ii})
    set(a.sub(ii), 'ylim', subYLim{ii}, 'ytick', subYTick{ii})
    if ii == 1; ylabel('voltage (µV)'); end
    
    plotTraces(elecResp, true, tracesMovie, cellInfo.gain, a.raw(ii))
    
    set(a.raw(ii), 'box', 'on', 'ylim', rawYLim{ii}, 'ytick', [], 'xtick', [])
    axes(a.raw(ii)); hold on
    tmpYLim = get(a.raw(ii), 'ylim');
    plot([1.85 1.85], tmpYLim(2)-0.1*(tmpYLim(2)-tmpYLim(1)) + [-100 0], 'k-', 'Linewidth', 2)
    plot([1.35 1.85], tmpYLim(2)-0.1*(tmpYLim(2)-tmpYLim(1))*[1 1], 'k-', 'Linewidth', 2)
    hold off
    
    clear datarun cellInfo shortName amp elecResp tracesMovie
end

%% plot color bar

colors = jet(101);

%plot vertical color scale bar

%figure('position', [100 100 100 500])
axes(a.colorBar)
hold on
nColors = size(colors, 1);
xPoly = [0.5 1 1 0.5];
for ii = 1:nColors
    yPoly = [(ii-0.5)/nColors (ii-0.5)/nColors (ii+0.5)/nColors (ii+0.5)/nColors];
    fill(xPoly, yPoly, colors(ii,:), 'EdgeColor', 'none')
end

colorBarLabels = 0:0.5:1;

for ii = 1:length(colorBarLabels)
    y = colorBarLabels(ii);
    %plot([0.4 0.5], [y y], 'k-')
    text(0.3, y, num2str(colorBarLabels(ii)), 'horizontalAlignment', 'right', 'verticalAlignment', 'middle')
end
set(gca, 'ylim', [-0.5 nColors+0.5]/nColors, 'xlim', [-1 1])
axis off
hold off


%% nested functions

    function plotRespCurve(elecResp, axesH, markedMovies, PW)
        data = zeros(3, length(elecResp.stimInfo.movieNos));
        data(1,:) = abs(elecResp.stimInfo.stimAmps);
        data(2,:) = elecResp.analysis.successRates;
        data(3,:) = elecResp.stimInfo.nPulses;
        starMov = false(1,size(data,2));
        for xx = length(elecResp.stimInfo.stimAmps): -1: 1
            if any(markedMovies == elecResp.stimInfo.movieNos(xx))
                starMov(xx) = true;
            end
            if isempty(elecResp.analysis.type{xx})
                data(:,xx) = [];
                starMov(xx) = [];
            end
        end
        
        %just in case erf fit isn't current
        elecResp = checkForUnfinishedAnalysis(elecResp, 100); %100 = bootstrap reps
        %[params, proj] = erfFitter(data, 2, -1);
        %thresh = -params(2)/params(1);
        

        
        xProj = 0:0.001:max(data(1,:))*2; %make sure this covers the whole displayed range
        proj = 0.5 + 0.5*erf(elecResp.analysis.erfParams(1)*xProj+elecResp.analysis.erfParams(2));

        %plot data and curve
        axes(axesH); hold on
        plot(xProj*PW, proj, 'k-')
        plot(data(1,~starMov)*PW, data(2,~starMov), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
        plot(data(1,find(starMov,1,'first'))*PW, data(2,find(starMov,1,'first')), 'o', 'MarkerFaceColor', [1 1 1], 'MarkerEdgeColor', 'k', 'MarkerSize', 5)
        plot(data(1,find(starMov,1,'last'))*PW, data(2,find(starMov,1,'last')), 'kx', 'MarkerSize', 8)
        hold off
        xlabel('charge amplitude (pC)')
    end

    function plotTraces(elecResp, raw, movieNo, gain, axesH)
        
        movInd = find(elecResp.stimInfo.movieNos == movieNo);
        
        dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
            elecResp.stimInfo.movieNos(movInd), 99999);
        
        
        centerChannel = elecResp.cells.recElec;
        %nPulses = elecResp.stimInfo.nPulses(movInd);
        latencies = [elecResp.analysis.latencies{movInd}];
                
        if raw == true
            dataToPlot = squeeze(dataTraces(:, centerChannel, :))*1000/gain;
        else %subtract estimated artifact
            subtractionVector = elecResp.analysis.estArtifact{movInd}(elecResp.cells.goodElecs == centerChannel, :)';
            elecResp.analysis.estArtifact{movInd}(elecResp.cells.goodElecs == elecResp.cells.recElec, :);
            subGrid = ndgrid(subtractionVector)';
            
            dataToPlot = (squeeze(dataTraces(:, centerChannel, :)) - subGrid(1:size(dataTraces,1),:))*1000/gain;
        end
        
        axes(axesH); hold on
        plot([0:0.05:2], dataToPlot(latencies==0, 2:42), 'k-') %plot sample 2 at t=0 because stim application occurs at sample 2
        plot([0:0.05:2], dataToPlot(latencies~=0, 2:42), 'r-')
        hold off
    end


end