function binMeans = lambdaVAxonPlotterBase(angles, lambdas, varargin)


p = inputParser;

p.addParamValue('makePlot', true, @islogical)
p.addParamValue('reflectData', false, @islogical)
p.addParamValue('binSize', 30, @isnumeric)
p.addParamValue('lamLimPlot', 1, @isnumeric)

p.parse(varargin{:})

params = p.Results;


if ~params.reflectData && max(angles)<=90
    binLims = 0:params.binSize:90;
else
    binLims = 0:params.binSize:180;
end

lam_binned = cell(1,length(binLims)-1);
binMeans = zeros(1, length(binLims)-1);

for ii = 1:length(binLims)-1
    lam_binned{ii} = lambdas(angles>binLims(ii) & angles<=binLims(ii+1));
    binMeans(ii) = mean(lam_binned{ii});
end


%% plotting

if params.makePlot
    figure; hold on    
    
    %plot grid for reference
    gridAmp = params.lamLimPlot;
    
    if params.reflectData
        plot(gridAmp*cosd(linspace(1,360,200)),      gridAmp*sind(linspace(1,360,200)),      '-', 'color', [0.8 0.8 0.8])
        plot(0.75*gridAmp*cosd(linspace(1,360,200)), 0.75*gridAmp*sind(linspace(1,360,200)), '-', 'color', [0.8 0.8 0.8])
        plot(0.5*gridAmp*cosd(linspace(1,360,200)),  0.5*gridAmp*sind(linspace(1,360,200)),  '-', 'color', [0.8 0.8 0.8])
        plot(0.25*gridAmp*cosd(linspace(1,360,200)), 0.25*gridAmp*sind(linspace(1,360,200)), '-', 'color', [0.8 0.8 0.8])
    elseif max(angles)<=90
        plot(gridAmp*cosd(linspace(1,90,100)),      gridAmp*sind(linspace(1,90,100)),      '-', 'color', [0.8 0.8 0.8])
        plot(0.75*gridAmp*cosd(linspace(1,90,100)), 0.75*gridAmp*sind(linspace(1,90,100)), '-', 'color', [0.8 0.8 0.8])
        plot(0.5*gridAmp*cosd(linspace(1,90,100)),  0.5*gridAmp*sind(linspace(1,90,100)),  '-', 'color', [0.8 0.8 0.8])
        plot(0.25*gridAmp*cosd(linspace(1,90,100)), 0.25*gridAmp*sind(linspace(1,90,100)), '-', 'color', [0.8 0.8 0.8])
    else
        plot(gridAmp*cosd(linspace(1,180,200)),      gridAmp*sind(linspace(1,180,200)),      '-', 'color', [0.8 0.8 0.8])
        plot(0.75*gridAmp*cosd(linspace(1,180,200)), 0.75*gridAmp*sind(linspace(1,180,200)), '-', 'color', [0.8 0.8 0.8])
        plot(0.5*gridAmp*cosd(linspace(1,180,200)),  0.5*gridAmp*sind(linspace(1,180,200)),  '-', 'color', [0.8 0.8 0.8])
        plot(0.25*gridAmp*cosd(linspace(1,180,200)), 0.25*gridAmp*sind(linspace(1,180,200)), '-', 'color', [0.8 0.8 0.8])
    end
    
    for ii = 1:length(binLims)
        plot([0 gridAmp*cosd(binLims(ii))], [0 gridAmp*sind(binLims(ii))], '-', 'color', [0.8 0.8 0.8])
        if params.reflectData
            plot([0 gridAmp*cosd(binLims(ii))], [0 -gridAmp*sind(binLims(ii))], '-', 'color', [0.8 0.8 0.8])
        end
    end
    
    %plot bins
    patchCoordsPos = zeros(2,1);
    patchCoordsNeg = zeros(2,1);
%     if binMeans(1)>0
%         plot([0 binMeans(1)], [0 0], 'k-')
%     else
%         plot([0 binMeans(1)], [0 0], 'r-')
%     end
    for ii = 1:length(binMeans)
        arc_th = linspace(binLims(ii), binLims(ii+1), 20);
        arc = abs(binMeans(ii))*[cosd(arc_th); sind(arc_th)];
        if binMeans(ii)>0
            patchCoordsPos = [patchCoordsPos arc]; %#ok<AGROW>
            patchCoordsNeg = [patchCoordsNeg [0;0]];
        else
            patchCoordsPos = [patchCoordsPos [0;0]]; %#ok<AGROW>
            patchCoordsNeg = [patchCoordsNeg arc];
        end
        %patchCoordsPos = [patchCoordsPos arc]; %#ok<AGROW>
    end
    if binLims(end) == 90
        patchCoordsPos = [patchCoordsPos [0; 0]];
        patchCoordsNeg = [patchCoordsNeg [0; 0]];
    end
    
    patch(patchCoordsPos(1,:), patchCoordsPos(2,:), [0.8 0.8 0.8])
    patch(patchCoordsNeg(1,:), patchCoordsNeg(2,:), [1 0.5 0.5])

    if params.reflectData
        patch(patchCoordsPos(1,:), -patchCoordsPos(2,:), [0.8 0.8 0.8])
        patch(patchCoordsNeg(1,:), -patchCoordsNeg(2,:), [1 0.5 0.5])
    end
    
    %plot data points
    for ii = 1:length(angles)
        if lambdas(ii)>0
            plot(lambdas(ii)*cosd(angles(ii)),     lambdas(ii)*sind(angles(ii)), 'k.')
            plot([0 lambdas(ii)*cosd(angles(ii))], [0 lambdas(ii)*sind(angles(ii))], 'k-')
        else
            plot(-lambdas(ii)*cosd(angles(ii)),     -lambdas(ii)*sind(angles(ii)), 'r.')
            plot([0 -lambdas(ii)*cosd(angles(ii))], [0 -lambdas(ii)*sind(angles(ii))], 'r-')
        end
        if params.reflectData
            if lambdas(ii)>0
                plot(lambdas(ii)*cosd(angles(ii)),     -lambdas(ii)*sind(angles(ii)), 'k.')
                plot([0 lambdas(ii)*cosd(angles(ii))], [0 -lambdas(ii)*sind(angles(ii))], 'k-')
            else
                plot(-lambdas(ii)*cosd(angles(ii)),     lambdas(ii)*sind(angles(ii)), 'r.')
                plot([0 -lambdas(ii)*cosd(angles(ii))], [0 lambdas(ii)*sind(angles(ii))], 'r-')
            end
        end
    end
    
    axis equal
    set(gca, 'xlim', [-gridAmp, gridAmp])
    if params.reflectData
        set(gca, 'ylim', [-gridAmp, gridAmp])
    else
        set(gca, 'ylim', [0, gridAmp])
    end
end

