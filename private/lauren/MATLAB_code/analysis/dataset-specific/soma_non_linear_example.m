clear all

errordlg('This script needs to be updated to work with the new erfFitter (now uses maximum likelihood fitting)')


pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data005';

neuronID = 800;
patternNos = [317 294]; %elec 49, elecs 46 + 49
nPatterns = length(patternNos);

cd(pathToData)
erfParams = cell(2);

colorScheme(1, :) = [27 117 187]/265; %blue
colorScheme(2, :) = [41 180 115]/265; %green
colorScheme(3, :) = [190 30 45]/265; %red

figure('position', [100 100 300 150])
hold on
for i = 1:nPatterns
    temp = load(['elecResp_n' num2str(neuronID) '_p' num2str(patternNos(i))]);
    elecResp = temp.elecResp;

    
    data = zeros(2, length(elecResp.stimInfo.movieNos));
    data(1,:) = abs(elecResp.stimInfo.stimAmps);
    data(2,:) = elecResp.analysis.successRates;
    for j = length(elecResp.stimInfo.stimAmps): -1: 1
        if isempty(elecResp.analysis.type{j})
            data(:,j) = [];
        end
    end

    erfParams = erfFitter(data, 2, -1, 'makePlot', 0);

    threshold = -erfParams(2)/erfParams(1);

    
    xProj = 0.1:0.001:1;
    projection = 0.5 + 0.5*erf(erfParams(1)*xProj+erfParams(2));
    
    current = plot(xProj, projection);
    set(findobj(current,'Type','line'),'Color', colorScheme(i,:), 'LineWidth', 2)
    xlabel('stimulus amplitude (uA)')
    ylabel('response probability')
end

projection = 0.5 + 0.5*erf(erfParams(1)*xProj+erfParams(2));
current = plot(xProj+1.5, projection, '--');
set(findobj(current,'Type','line'),'Color', colorScheme(3,:), 'LineWidth', 2)
set(gca, 'xlim', [0.1 2.5])


hold off

figure
plotEi61('/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data004-NW/data004-NW.ei', 800, gca, [49 46])