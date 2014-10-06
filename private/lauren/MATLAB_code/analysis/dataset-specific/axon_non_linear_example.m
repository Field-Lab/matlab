clear all

errordlg('This script needs to be updated to work with the new erfFitter (now uses maximum likelihood fitting)')


pathToData = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data010';

neuronID = 559; 

movieNos = 15:26;
patternNo = 49; %just to get stim amps

elec35Resp = [0 0 0 0.89 0.99 0.99 1 1 1 1 1 1];
elec37Resp = [0 0 0 0 0 0 0 0 0 0.02 1 0.99];
bothElecsResp = [0 0 0 0.01 0.1 0.96 1 1 1 1 1 1];

responses{1} = elec35Resp;
responses{2} = elec37Resp;
responses{3} = bothElecsResp;

%assumes that stim amps are same for 2008-08-26-0/data010 and 2008-08-27-4/data006
stimAmps = zeros(1, length(movieNos));
for i = 1:length(movieNos)
    stimAmps(i) = getStimAmps(pathToData, patternNo, movieNos(i));
end


cd(pathToData)
erfParams = cell(3);

colorScheme(1, :) = [27 117 187]/265; %blue
colorScheme(2, :) = [190 30 45]/265; %red
colorScheme(3, :) = [41 180 115]/265; %green

figure('position', [100 100 200 150])
hold on
for i = 1:3
    disp(num2str(i))

    data = [abs(stimAmps); responses{i}];
    erfParams = erfFitter(data, 2, 1, 'makePlot', 0);
    threshold = -erfParams(2)/erfParams(1);
    
    

    xProj = 0.4:0.001:1;
    projection = 0.5 + 0.5*erf(erfParams(1)*xProj+erfParams(2));
    
    current = plot(xProj, projection);
    set(findobj(current,'Type','line'),'Color', colorScheme(i,:), 'LineWidth', 2)
    xlabel('stimulus amplitude (uA)')
    ylabel('response probability')
end
hold off

%figure
%plotEi61('/snle/lab/Experiments/Array/Analysis/2008-08-27-4/data006/data006.ei', 227, gca, [35 37])