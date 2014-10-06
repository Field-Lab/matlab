function COM = getCOM(ei, amps, threshold, makePlot, checkForAxons, id, varargin)

% computes 'center of mass' of ei, based on ei values above the specified threshold
%
%
%

p = inputParser;

p.addParamValue('plotAxes', [], @ishandle)
p.addParamValue('plotFigAxonCheck', [], @ishandle)
p.addParamValue('axonElecs', [], @isnumeric)

p.parse(varargin{:})

params = p.Results;

%%

[xCoords yCoords] = getElectrodeCoords61();


threshAmp = max(amps)*threshold;

COM = zeros(2,1);
for i = 1:length(amps)
    if amps(i) < threshAmp || any(params.axonElecs == i)
        amps(i) = 0;
    else
        COM(1) = COM(1) + amps(i)*xCoords(i);
        COM(2) = COM(2) + amps(i)*yCoords(i);
    end
end

COM = COM/sum(amps);

if makePlot
    if isempty(params.plotAxes)
        figure
    else
        axes(params.plotAxes); cla(params.plotAxes)
    end
    hold on
    for i = 1:length(amps)
        if amps(i)>0
            plot(xCoords(i), yCoords(i), 'ok','MarkerSize', ceil(amps(i)*20), 'MarkerFaceColor', 'k')
        end
    end
    plot(COM(1), COM(2), 'or', 'MarkerSize', 10)
    
    %plot array outline
    plot([0 8.6603 8.6603 0 -8.6603 -8.6603 0], [10 5 -5 -10 -5 5 10], 'k-')
    hold off
    
    set(gca, 'XLim', [-10 10], 'YLim', [-11 11])
    axis equal
    axis off
    if exist('id', 'var')
        title(['neuron ' num2str(id)])
    end
end




%check if any additional axonal signals need to be excluded
if checkForAxons
    arrayWidth = max(xCoords)*2*1.2;
    arrayHeight = max(yCoords)*2*1.2;
    if isempty(params.plotAxes)
        figure('position', [100 100 400 350])
    else
        figure(params.plotFigAxonCheck); clf(params.plotFigAxonCheck)
    end
    
    if exist('id', 'var')
        title(['neuron ' num2str(id)])
        axis off
    end
    
    for i = 1:64
        if ~(i==9||i==25||i==57)
            axes('position', [xCoords(i)/arrayWidth + 0.45, yCoords(i)/arrayHeight + 0.45, 0.1, 0.1])
            plot(squeeze(ei(1,i+1,:)), 'k')            
            if amps(i) > 0
                set(gca, 'box', 'on')
            elseif any(params.axonElecs == i)
                plot(squeeze(ei(1,i+1,:)), 'r')
            else
                axis off
            end
            set(gca, 'ylim', [-1 1], 'xlim', [1 60], 'xtick', [], 'ytick', [])
        end
    end
end