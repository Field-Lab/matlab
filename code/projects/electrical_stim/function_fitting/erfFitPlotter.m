function erfFitPlotter(elecResps, axesHandle, plotColors, logScale)

% makes plot of multiple response curves on one graph, with best-fit error function, in specified
% colors
%
% for use with axonStimFigurePlotter and somaStimFigurePlotter (someStimPlotter needs to updated
% first)

if ~exist('logScale', 'var')
    logScale = 0;
end

nCurves = length(elecResps);

%axes(axesHandle)

hold on
for i = 1:nCurves
    
    data = zeros(2, length(elecResps{i}.stimInfo.stimAmps));
    data(1,:) = abs(elecResps{i}.stimInfo.stimAmps);
    data(2,:) = elecResps{i}.analysis.successRates;
    for j = length(elecResps{i}.stimInfo.stimAmps): -1: 1
        if isempty(elecResps{i}.analysis.type{j})
            data(:,j) = [];
        end
    end

    if logScale
        semilogx(data(1,:), data(2,:),'o','MarkerEdgeColor', plotColors(i,:),'MarkerFaceColor', plotColors(i,:))
    else
        plot(data(1,:), data(2,:),'o','MarkerEdgeColor', plotColors(i,:),'MarkerFaceColor', plotColors(i,:))
    end

    xProj = min(abs(data(1,:))) : 0.01 : max(abs(data(1,:)));
    projection = 0.5 + 0.5*erf(elecResps{i}.analysis.erfParams(1)*xProj + elecResps{i}.analysis.erfParams(2));
    
    if logScale
        plot(xProj, projection, 'Color', plotColors(i,:));
    else
        plot(xProj, projection, 'Color', plotColors(i,:));
    end
end
hold off