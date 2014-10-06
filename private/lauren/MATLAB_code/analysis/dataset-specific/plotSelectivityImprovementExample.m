function plotSelectivityImprovementExample(CI, bestPInd)


nPatterns = length(CI(1).patternNos);

meanRelAmps = cellfun(@(x) mean(x,2), CI(1).details.actualRelAmps, 'uniformOutput', false);
meanRelAmps = cell2mat(meanRelAmps')';

meanPAmps = cellfun(@(x) mean(x,2), CI(1).details.pAmps, 'uniformOutput', false);
meanPAmps = cell2mat(meanPAmps')';

meanSAmps = cellfun(@(x) mean(x,2), CI(1).details.sAmps, 'uniformOutput', false);
meanSAmps = cell2mat(meanSAmps')';

pAloneBin = CI(1).details.pAloneBin & meanPAmps < 0;
sAloneBin = CI(1).details.sAloneBin & all(meanSAmps <= 0, 2);

valid = cell(2,1);
valid{1} = true(nPatterns, 1); valid{2} = true(nPatterns, 1);
if isfield((CI(1).details),'unfitTrip') %LG mod
    valid{1}(CI(1).details.unfitTrip | CI(1).details.unfitPairs | CI(1).details.excludeFromFits |...
        CI(1).details.pAloneBin | CI(1).details.sAloneBin | sum(meanRelAmps~=0, 2)>2) = false;
    valid{2}(CI(2).details.unfitTrip | CI(2).details.unfitPairs | CI(2).details.excludeFromFits |...
        CI(1).details.pAloneBin | CI(1).details.sAloneBin | sum(meanRelAmps~=0, 2)>2) = false;
else
    valid{1}( CI(1).details.unfitPairs | CI(1).details.excludeFromFits |...
        CI(1).details.pAloneBin | CI(1).details.sAloneBin | sum(meanRelAmps~=0, 2)>2) = false;
    valid{2}( CI(2).details.unfitPairs | CI(2).details.excludeFromFits |...
        CI(1).details.pAloneBin | CI(1).details.sAloneBin | sum(meanRelAmps~=0, 2)>2) = false;
end
%exclude patterns outside primary positive region
for ii = 1:nPatterns
    for jj = 1:2
        for kk = 1:6
            relAmp = meanRelAmps(ii,kk);
            mBound = CI(jj).details.mValidBound(:,kk);
            if (mBound(1) && relAmp > 0 && (1/relAmp) < mBound(1)) || (mBound(2) && relAmp < 0 && abs(1/relAmp) < abs(mBound(2)))
                valid{jj}(ii) = false;
            end
        end
    end
end


predThresh = zeros(nPatterns, 2);

%estimated threshold for each valid combination, based on full model
%mean relative amplitudes on component pairs
for ii = 1:nPatterns
    if valid{1}(ii) && valid{2}(ii)
        sInds = find(meanRelAmps(ii,:));
        for jj = 1:2 %loop through cells
            yInt = CI(jj).details.tp; %model primary-alone threshold
            shiftSum = 0; %predicted shift from primary-alone threshold
            for kk = 1:length(sInds)
                relAmp = meanRelAmps(ii,sInds(kk));
                
                %find I_i at intercept between fit line of model (yIntercept-lambda_i*I_i) and line representing current combinations
                %used in the pattern (1/relAmps)*I_i (this intercept represents the threshold predicted by the model)
                lam = CI(jj).details.lambdas(sInds(kk));
                I_i = yInt/((1/relAmp)+lam);
                
                %shift = -1*lambda_i*I_i
                shiftSum = shiftSum - lam*I_i;
            end
            predThresh(ii,jj) = yInt + shiftSum;
        end
    end
end

threshDiffPred = predThresh(:,1)-predThresh(:,2);
threshDiffData = CI(1).thresholds - CI(2).thresholds;

threshDiffPAloneData = threshDiffData(pAloneBin);
threshDiffPAlonePred = CI(1).details.tp - CI(2).details.tp;

threshDiffData(~(valid{1} & valid{2})) = -inf;
threshDiffPred(~(valid{1} & valid{2})) = -inf;

[~, top10Pred] = sort(threshDiffPred, 1, 'descend');
[~, top10Data] = sort(threshDiffData, 1, 'descend');
top10Pred = top10Pred(1:10);
top10Data = top10Data(1:10);

threshDiffData(~(valid{1} & valid{2})) = [];
threshDiffPred(~(valid{1} & valid{2})) = [];

[~, predOrd] = sort(threshDiffPred, 1, 'descend');

figure('position', [100 100 500 180]); hold on
plot(threshDiffPred(predOrd), 'k.')
plot(threshDiffData(predOrd), 'r.')
for ii = 1:length(predOrd)
    plot(ii*[1 1], [threshDiffPred(predOrd(ii)) threshDiffData(predOrd(ii))], 'r-')
end
plot(-2, threshDiffPAlonePred, 'k.')
plot(-2, threshDiffPAloneData, 'b.')
plot(-2*[1 1], [threshDiffPAlonePred threshDiffPAloneData], 'b-')
ylabel('threshold difference (µA)')
xlabel('spatial pattern')
set(gca, 'xlim', [-3 length(predOrd)+1])


%% plot response curves for top 3 patterns 

[xCoords yCoords] = getElectrodeCoords61();
%center on primary electrode
xCoords = xCoords - xCoords(CI(jj).pElec);
yCoords = yCoords - yCoords(CI(jj).pElec);
sElecCoords(1,:) = xCoords(CI(jj).details.sElecs);
sElecCoords(2,:) = yCoords(CI(jj).details.sElecs);

cellColors(1,:) = 0*[0.25 0.25 0.25];
cellColors(2,:) = [90 156 0]/255;

% best according to data
elecRespData = cell(2,3);
figure('position', [100 100 1000 500])

erfParamsAll = cell(2,3);

for ii = 1:3
    axes('units', 'pixels', 'position', [50+300*(ii-1) 20 125 90])
    hold on
    mu = zeros(1,2); sig = zeros(1,2);
    for jj = 1:2
        tmp = load([CI(jj).pathToData filesep 'elecResp_n' num2str(CI(jj).id) '_p' num2str(CI(jj).patternNos(top10Data(ii)))]);
        elecRespData{jj,ii} = tmp.elecResp;
        
        data = zeros(3, length(elecRespData{jj,ii}.stimInfo.movieNos));
        data(1,:) = abs(elecRespData{jj,ii}.stimInfo.stimAmps);
        data(2,:) = elecRespData{jj,ii}.analysis.successRates;
        data(3,:) = elecRespData{jj,ii}.stimInfo.nPulses;
        
        valData = ~cellfun(@isempty, elecRespData{jj,ii}.analysis.type);
        data(:,~valData) = [];
        
        erfParams = elecRespData{jj,ii}.analysis.erfParams;
        xProj = 0:0.01:0.8;
        projection = 0.5 + 0.5*erf(erfParams(1)*xProj + erfParams(2));
        
        plot(data(1,:), data(2,:), '.', 'markerEdgeColor', cellColors(jj,:))
        plot(xProj, projection, '-', 'color', cellColors(jj,:))
        
        mu(jj) = -erfParams(2)/erfParams(1);
        sig(jj) = 1/(sqrt(2)*erfParams(1));
    end
    
    %amplitude required to achieve 0.9 probability for cell 2
    amp90 = norminv(0.9, mu(2), sig(2));
    
    %cell 1 response probability at this amplitude:
    cell1RespProb = normcdf(amp90, mu(1), sig(1));
    
    
    plot([0 0.8], [0.9 0.9])
    title(['pattern ' num2str(CI(jj).patternNos(top10Data(ii))), 10, 'cell 1 prob = ' num2str(cell1RespProb)])

    set(gca, 'xlim', [0 0.8], 'xtick', [0 0.4 0.8], 'ytick', [0 0.5 1])
    
    %plot stimulus
    ra = meanRelAmps(top10Data(ii),:);
    axes('units', 'pixels', 'position', [50+300*(ii-1) 350 50 50])
    hold on
    plot(0, 0, 'k.')
    plot([0 0], [0 1], 'k-')
    plot(sElecCoords(1,:), sElecCoords(2,:), 'k.')
    for jj = 1:6
        if ra(jj)>0
            plot(sElecCoords(1,jj)*[1 1], sElecCoords(2,jj)+[0 abs(ra(jj))], 'k-')
        elseif ra(jj)<0
            plot(sElecCoords(1,jj)*[1 1], sElecCoords(2,jj)+[0 abs(ra(jj))], 'r-')
        end
    end
    set(gca, 'xlim', [-5 5], 'ylim', [-5 5])
end


% best according to model
elecRespPred = cell(2,3); %cell x pattern
figure('position', [100 100 1000 500])

for ii = 1:3
    axes('units', 'pixels', 'position', [50+300*(ii-1) 20 125 90])
    hold on
    mu = zeros(1,2); sig = zeros(1,2);
    for jj = 1:2
        tmp = load([CI(jj).pathToData filesep 'elecResp_n' num2str(CI(jj).id) '_p' num2str(CI(jj).patternNos(top10Pred(ii)))]);
        elecRespPred{jj,ii} = tmp.elecResp;
        
        data = zeros(3, length(elecRespPred{jj,ii}.stimInfo.movieNos));
        data(1,:) = abs(elecRespPred{jj,ii}.stimInfo.stimAmps);
        data(2,:) = elecRespPred{jj,ii}.analysis.successRates;
        data(3,:) = elecRespPred{jj,ii}.stimInfo.nPulses;
        
        valData = ~cellfun(@isempty, elecRespPred{jj,ii}.analysis.type);
        data(:,~valData) = [];
        
        erfParams = elecRespPred{jj,ii}.analysis.erfParams;
        xProj = 0:0.01:0.8;
        projection = 0.5 + 0.5*erf(erfParams(1)*xProj + erfParams(2));
        
        plot(data(1,:), data(2,:), '.', 'markerEdgeColor', cellColors(jj,:))
        plot(xProj, projection, '-', 'color', cellColors(jj,:))
        
        %convert parameters to mean and stdev of cumulative gaussian
        mu(jj) = -erfParams(2)/erfParams(1);
        sig(jj) = 1/(sqrt(2)*erfParams(1));
    end
    
    %amplitude required to achieve 0.9 probability for cell 2
    amp90 = norminv(0.9, mu(2), sig(2));
    
    %cell 1 response probability at this amplitude:
    cell1RespProb = normcdf(amp90, mu(1), sig(1));
    
    plot([0 0.8], [0.9 0.9])
    title(['pattern ' num2str(CI(jj).patternNos(top10Pred(ii))), 10, 'cell 1 prob = ' num2str(cell1RespProb)])
    set(gca, 'xlim', [0 0.8], 'xtick', [0 0.4 0.8], 'ytick', [0 0.5 1])
    
    %plot stimulus
    ra = meanRelAmps(top10Pred(ii),:);
    axes('units', 'pixels', 'position', [50+300*(ii-1) 350 50 50])
    hold on
    plot(0, 0, 'k.')
    plot([0 0], [0 1], 'k-')
    plot(sElecCoords(1,:), sElecCoords(2,:), 'k.')
    for jj = 1:6
        if ra(jj)>0
            plot(sElecCoords(1,jj)*[1 1], sElecCoords(2,jj)+[0 abs(ra(jj))], 'k-')
        elseif ra(jj)<0
            plot(sElecCoords(1,jj)*[1 1], sElecCoords(2,jj)+[0 abs(ra(jj))], 'r-')
        end
    end
    set(gca, 'xlim', [-5 5], 'ylim', [-5 5])
end


%primary-alone
figure('position', [100 100 350 500])
axes('units', 'pixels', 'position', [50 20 125 90])
hold on
mu = zeros(1,2); sig = zeros(1,2);
for jj = 1:2
    tmp = load([CI(jj).pathToData filesep 'elecResp_n' num2str(CI(jj).id) '_p' num2str(CI(jj).patternNos(pAloneBin))]);
    er = tmp.elecResp;
    
    data = zeros(3, length(er.stimInfo.movieNos));
    data(1,:) = abs(er.stimInfo.stimAmps);
    data(2,:) = er.analysis.successRates;
    data(3,:) = er.stimInfo.nPulses;
    
    valData = ~cellfun(@isempty, er.analysis.type);
    data(:,~valData) = [];
    
    erfParams = er.analysis.erfParams;
    xProj = 0:0.01:0.8;
    projection = 0.5 + 0.5*erf(erfParams(1)*xProj + erfParams(2));
    
    plot(data(1,:), data(2,:), '.', 'markerEdgeColor', cellColors(jj,:))
    plot(xProj, projection, '-', 'color', cellColors(jj,:))
    
    mu(jj) = -erfParams(2)/erfParams(1);
    sig(jj) = 1/(sqrt(2)*erfParams(1));
end

%amplitude required to achieve 0.9 probability for cell 2
amp90 = norminv(0.9, mu(2), sig(2));

%cell 1 response probability at this amplitude:
cell1RespProb = normcdf(amp90, mu(1), sig(1));

plot([0 0.8], [0.9 0.9])
title(['pattern ' num2str(CI(jj).patternNos(pAloneBin)), 10, 'cell 1 prob = ' num2str(cell1RespProb)])
set(gca, 'xlim', [0 0.8], 'xtick', [0 0.4 0.8], 'ytick', [0 0.5 1])

%plot stimulus
axes('units', 'pixels', 'position', [50 350 50 50])
hold on
plot(0, 0, 'k.')
plot(sElecCoords(1,:), sElecCoords(2,:), 'k.')
plot([0 0], [0 1], 'k-')
set(gca, 'xlim', [-5 5], 'ylim', [-5 5])

%secondary-alone
sAloneInd = find(sAloneBin);
for ii = 1:6
    figure('position', [100 100 350 500])
    axes('units', 'pixels', 'position', [50 20 125 90])
    hold on
    for jj = 1:2
        tmp = load([CI(jj).pathToData filesep 'elecResp_n' num2str(CI(jj).id) '_p' num2str(CI(jj).patternNos(sAloneInd(ii)))]);
        er = tmp.elecResp;
        
        data = zeros(3, length(er.stimInfo.movieNos));
        data(1,:) = abs(er.stimInfo.stimAmps);
        data(2,:) = er.analysis.successRates;
        data(3,:) = er.stimInfo.nPulses;
        
        valData = ~cellfun(@isempty, er.analysis.type);
        data(:,~valData) = [];
        
        erfParams = er.analysis.erfParams;
        xProj = 0:0.01:1.6;
        projection = 0.5 + 0.5*erf(erfParams(1)*xProj + erfParams(2));
        
        plot(data(1,:), data(2,:), '.', 'markerEdgeColor', cellColors(jj,:))
        plot(xProj, projection, '-', 'color', cellColors(jj,:))
    end
    plot([0 1.6], [0.9 0.9])
    title(['pattern ' num2str(CI(jj).patternNos(sAloneInd(ii)))])
    set(gca, 'xlim', [0 1.6], 'xtick', [0 0.4 0.8 1.2 1.6], 'ytick', [0 0.5 1])
    
    %plot stimulus
    sInd = find(isinf(meanRelAmps(sAloneInd(ii),:)));
    axes('units', 'pixels', 'position', [50 350 50 50])
    hold on
    plot(0, 0, 'k.')
    plot(sElecCoords(1,:), sElecCoords(2,:), 'k.')
    plot(sElecCoords(1,sInd)*[1 1], sElecCoords(2,sInd)+[0 1], 'k-')
    set(gca, 'xlim', [-5 5], 'ylim', [-5 5])
end
keyboard


%%
        
%outsideReg{1}

colValid = [0 0 0; 0 0 1];
colInvalid = [0.5 0.5 0.5;  0.5 0.5 1];

% plot data and planar fit for each pair of secondary electrodes

x = cell(2,1); y = cell(2,1); z = cell(2,1); lam = cell(2,1);

figure('position', [100 100 900 600])
hAxes{1} = axes('position', [0.05 0.55 0.25 0.4]);
hAxes{2} = axes('position', [0.35 0.55 0.25 0.4]);
hAxes{3} = axes('position', [0.65 0.55 0.25 0.4]);
hAxes{4} = axes('position', [0.05 0.05 0.25 0.4]);
hAxes{5} = axes('position', [0.35 0.05 0.25 0.4]);
hAxes{6} = axes('position', [0.65 0.05 0.25 0.4]);


for ii = 1:6
    x{1} = zeros(nPatterns,1); x{2} = zeros(nPatterns,1);
    y{1} = zeros(nPatterns,1); y{2} = zeros(nPatterns,1);
    z{1} = zeros(nPatterns,1); z{2} = zeros(nPatterns,1);
    
    if ii~=6
        sElecPair = [ii ii+1];
    else
        sElecPair = [6 1];
    end
        
    ra1 = zeros(nPatterns, 1);
    ra2 = zeros(nPatterns, 2);
    for jj = 1:nPatterns
        ra1(jj) = meanRelAmps(jj, sElecPair(1));
        ra2(jj) = meanRelAmps(jj, sElecPair(2));
        raRest = meanRelAmps(jj, setdiff(1:6, sElecPair));
        
        if ((ra1(jj) || ra2(jj)) && ~any(raRest) && ~CI(1).details.sAloneBin(jj)) ||...
                (CI(1).details.pAloneBin(jj) && CI(1).details.pAmps{jj}(end)<0) %at least one of the sElecs is involved or cathodal primary-alone
            for kk = 1:2 %loop through cells
                x{kk}(jj) = CI(kk).thresholds(jj)*ra1(jj);
                y{kk}(jj) = CI(kk).thresholds(jj)*ra2(jj);
                z{kk}(jj) = CI(kk).thresholds(jj);
            end
        end
    end
    
    axes(hAxes{ii}); hold on
    for jj = 1:nPatterns
        if valid{1}(jj) && valid{2}(jj) && z{1}(jj) ~= 0
            plot3([x{1}(jj) x{2}(jj)], [y{1}(jj) y{2}(jj)], [z{1}(jj) z{2}(jj)], '-', 'color', [0 0 0])
            if any(bestPInd == jj)
                plot3([x{1}(jj) x{2}(jj)], [y{1}(jj) y{2}(jj)], [z{1}(jj) z{2}(jj)], '-', 'color', [1 0 0], 'linewidth', 2)
            end
        end
    end
    for kk = 1:2
        toPlotValid = valid{kk} & z{kk}~=0;
        toPlotInvalid = ~valid{kk} & z{kk}~=0;

        plot3(x{kk}(toPlotValid), y{kk}(toPlotValid), z{kk}(toPlotValid), '.',...
            'markerEdgeColor', colValid(kk,:), 'markerFaceColor', colValid(kk,:))
        plot3(x{kk}(toPlotInvalid), y{kk}(toPlotInvalid), z{kk}(toPlotInvalid), '.',...
            'markerEdgeColor', colInvalid(kk,:), 'markerFaceColor', colInvalid(kk,:))
    end
    
    %find the optimal pattern from this electrode pair: should be one
    %corner (within valid region)
    lam1Cell1 = CI(1).details.lambdas(sElecPair(1));
    lam2Cell1 = CI(1).details.lambdas(sElecPair(2));
    
    lam1Cell2 = CI(2).details.lambdas(sElecPair(1));
    lam2Cell2 = CI(2).details.lambdas(sElecPair(2));

    relAmpMaxVal1 = max([CI(1).details.mValidBound(1, sElecPair(1)) CI(2).details.mValidBound(1, sElecPair(1))]);
    relAmpMinVal1 = -1*max(abs([CI(1).details.mValidBound(2, sElecPair(1)) CI(2).details.mValidBound(2, sElecPair(1))]));
    
    relAmpMaxVal2 = max([CI(1).details.mValidBound(1, sElecPair(2)) CI(2).details.mValidBound(1, sElecPair(2))]);
    relAmpMinVal2 = -1*max(abs([CI(1).details.mValidBound(2, sElecPair(2)) CI(2).details.mValidBound(2, sElecPair(2))]));

    
%     if lam1Cell1 > lam1Cell2
%         optDir1 = 'pos';
%         optRelAmp1 = max()
%     else
%         optDir1 = 'neg';
%     end
%     
%     if lam2Cell1 > lam2Cell2
%         optDir2 = 'pos';
%     else
%         optDir2 = 'neg';
%     end
%         
    
    
    %plot valid region boundaries
    xRange = [min([min(x{1}), min(x{2})]) max([max(x{1}), max(x{2})])];
    yRange = [min([min(y{1}), min(y{2})]) max([max(y{1}), max(y{2})])];
    
    mx{1} = CI(1).details.mValidBound(:, sElecPair(1));
    mx{2} = CI(2).details.mValidBound(:, sElecPair(1));
    
    my{1} = CI(1).details.mValidBound(:, sElecPair(2));
    my{2} = CI(2).details.mValidBound(:, sElecPair(2));
    
    gridInd = 0.2;
    
    for kk = 1:2
        if mx{kk}(1) ~= 0
            [xMesh yMesh] = meshgrid(0:gridInd:xRange(2), yRange(1):gridInd:yRange(2));
            mesh(xMesh, yMesh, xMesh*mx{kk}(1), 'faceColor', 'none', 'edgeColor', colValid(kk,:)) %cathodal (positive on axis)
        end
        
        if mx{kk}(2) ~= 0
            [xMesh yMesh] = meshgrid(xRange(1):gridInd:0, yRange(1):gridInd:yRange(2));
            mesh(xMesh, yMesh, xMesh*mx{kk}(2), 'faceColor', 'none', 'edgeColor', colValid(kk,:)) %anodal (negative on axis)
        end
        
        if my{kk}(1) ~= 0
            [xMesh yMesh] = meshgrid(xRange(1):gridInd:xRange(2), 0:gridInd:yRange(2));
            mesh(xMesh, yMesh, yMesh*my{kk}(1), 'faceColor', 'none', 'edgeColor', colValid(kk,:)) %cathodal (positive on axis)
        end
        
        if my{kk}(2) ~= 0
            [xMesh yMesh] = meshgrid(xRange(1):gridInd:xRange(2), yRange(1):gridInd:0);
            mesh(xMesh, yMesh, yMesh*my{kk}(2), 'faceColor', 'none', 'edgeColor', colValid(kk,:)) %anodal (negative on axis)
        end
    end
        
%     %valid region boundaries for cell 2
%     [xMesh yMesh] = meshgrid(0:gridInd:xRange(2), yRange(1):gridInd:yRange(2));
%     mesh(xMesh, yMesh, xMesh*mx2(1), 'faceColor', 'none', 'edgeColor', colValid(2,:)) %cathodal (positive on axis)
% 
%     [xMesh yMesh] = meshgrid(xRange(1):gridInd:0, yRange(1):gridInd:yRange(2));
%     mesh(xMesh, yMesh, xMesh*mx2(2), 'faceColor', 'none', 'edgeColor', colValid(2,:)) %anodal (negative on axis)
%     
%     [xMesh yMesh] = meshgrid(xRange(1):gridInd:xRange(2), 0:gridInd:yRange(2));
%     mesh(xMesh, yMesh, yMesh*my2(1), 'faceColor', 'none', 'edgeColor', colValid(2,:)) %cathodal (positive on axis)
% 
%     [xMesh yMesh] = meshgrid(xRange(1):gridInd:xRange(2), yRange(1):gridInd:0);
%     mesh(xMesh, yMesh, yMesh*my2(2), 'faceColor', 'none', 'edgeColor', colValid(2,:)) %anodal (negative on axis)

    
    set(hAxes{ii}, 'xlim', [-1 1], 'ylim', [-1 1], 'zlim', [0 1])
    axis equal
end


