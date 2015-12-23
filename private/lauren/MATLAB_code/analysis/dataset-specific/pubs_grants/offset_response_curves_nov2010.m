%% plots multiple response curves after normalizing x-axis to threshold and offsetting curves for
% clarity
%
% used in grant resubmission Nov. 2010
%
%

cd /snle/lab/Experiments/Array/Analysis/2008-08-27-2/data002/

clear all

recalcAll = 0; %redo all response curve fits and standard deviation calculations
keepLogBased = 0; %calculates log-based erf fit, threshold, and thresh standard deviations in addition to linear-based
nBootstrapReps = 100;


%% loads up data


midgetFilenames = {'n18_p2', 'n91_p7', 'n92_p6', 'n167_p12', 'n182_p13', 'n256_p20', 'n274_p19',...
    'n332_p23', 'n348_p24', 'n379_p26', 'n395_p27', 'n512_p35', 'n541_p37', 'n616_p42', 'n646_p44',...
    'n676_p46', 'n677_p46', 'n766_p53', 'n782_p54', 'n811_p55', 'n858_p58', 'n872_p61', 'n886_p60',...
    'n931_p63'};


nMidgets = length(midgetFilenames);


for i = 1:nMidgets
    midgetFilenames{i} = ['elecResp_' midgetFilenames{i}];
end

midgetData =        cell(nMidgets, 1); %stores responses and stimulus amplitudes
midgetParams      = cell(nMidgets, 1);
midgetSlopes = zeros(nMidgets, 1);


for i = 1:nMidgets
    temp = load(midgetFilenames{i});
    elecResp = temp.elecResp;
    
    elecResp = checkForUnfinishedAnalysis(elecResp, nBootstrapReps, 'recalcAll', recalcAll, 'keepLogBased', keepLogBased);
    save(midgetFilenames{i}, 'elecResp')
    
    nMovies = length(elecResp.stimInfo.movieNos);
    midgetData{i} = zeros(2, nMovies);
    midgetData{i}(1,:) = abs(elecResp.stimInfo.stimAmps);
    midgetData{i}(2,:) = elecResp.analysis.successRates;
    %removes stimulus amplitudes that have not been analyzed
    for j = length(elecResp.stimInfo.stimAmps): -1: 1
        if isempty(elecResp.analysis.type{j})
            midgetData{i}(:,j) = [];
        end
    end
    midgetParams{i} = elecResp.analysis.erfParams;
    midgetSlopes(i) = -sqrt(2)*elecResp.analysis.erfParams(2);
    
    midgetThresholds(i) = elecResp.analysis.threshold;
end

%% makes plot

figure('position', [100 100 325 200], 'color', [1 1 1])
hold on

if 1
    n = 0;
    for i = 1:nMidgets
        if max(midgetData{i}(2,:)) > 0.95
            n = n+1;

            xProj = midgetData{i}(1,1):0.01:midgetData{i}(1,end);
            projection = 0.5 + 0.5*erf(midgetParams{i}(1)*xProj + midgetParams{i}(2));

            plot(midgetData{i}(1,:)/midgetThresholds(i) + n*0.6, midgetData{i}(2,:), 'k.') %used in grant
            plot(xProj/midgetThresholds(i) + n*0.6, projection, 'k-'); %used in grant
            
%             for jj = 1:size(midgetData{i},2)
%                 xVals = (midgetData{i}(1,jj)/midgetThresholds(i) + n*0.6)*ones(1,2);
%                 yVals = 100*[midgetData{i}(2,jj), 0.5+0.5*erf(midgetParams{i}(1)*midgetData{i}(1,jj) + midgetParams{i}(2))];
%                 plot(xVals, yVals, 'k-')
%             end
            
            
            %    plot(midgetData{i}(1,:)/midgetThresholds(i), midgetData{i}(2,:)*100, '.', 'MarkerEdgeColor', midgetColor,...
            %        'MarkerFaceColor', midgetColor, 'MarkerSize', 10)
            %plot(midgetData{i}(1,:)/midgetThresholds(i), midgetData{i}(2,:)*100, '-', 'color', midgetColor)
            %plot(xProj/midgetThresholds(i), projection*100, 'Color', midgetColor);
        end
    end
    set(gca, 'XLim', [0.5 7])
end


