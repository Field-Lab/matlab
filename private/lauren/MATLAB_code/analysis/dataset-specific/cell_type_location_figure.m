%% get cell type properties from all cells in each class in data004 of this piece (time courses and
% average RF diameters)


% gets ratio of areas in transformation from monitor space to array space (from a dataset with the
% same set-up)
datarun = load_data('2010-03-05-2/data000/data000');            %      519, 30 µm, 6.5x from below

% load datarun
datarun = load_index(datarun);
datarun = load_sta(datarun,'load_sta',[]);
datarun = load_params(datarun);
datarun = load_ei(datarun,[]);
datarun = get_sta_fits_from_vision(datarun);

% load alignment information from disk
datarun = load_monitor_alignment(datarun);

arrayAreaMonitorSpace = polyarea(datarun.piece.corners(1:end-1, 1), datarun.piece.corners(1:end-1, 2));

arrayCornerPos = datarun.ei.position(datarun.piece.corner_electrodes, :);
arrayAreaElectrodeSpace = polyarea(arrayCornerPos(:,1), arrayCornerPos(:,2));

areaRatio = arrayAreaElectrodeSpace/arrayAreaMonitorSpace;


%% loads dataset with same cells used in center of mass plot, but with better time courses and rf
% fits

datarun = load_data('2008-08-27-5/data004/data004/data004');
datarun = load_index(datarun);

datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun, 'save_rf', 1);
datarun = load_ei(datarun, []);
datarun = load_neurons(datarun);
datarun = get_sta_summaries(datarun, 'all');

datarun = get_sta_fits_from_vision(datarun,'all');

% for ii = 1:5
%     figure
%     plot_rf_summaries(datarun, {ii})
% end

cellTypeIndices = cell(5,1);
for ii = 1:5
    cellTypeIndices{ii} = get_cell_indices(datarun, {ii});
end


%removes crappy cells from on midget mosaic
for ii = length(cellTypeIndices{3}): -1 :1
    if any([2016 5210 5780 7373] == datarun.cell_ids(cellTypeIndices{3}(ii)))
        cellTypeIndices{3}(ii) = [];
    end
end


for ii = 1:5
    figure
    plot_rf_summaries(datarun, datarun.cell_ids(cellTypeIndices{ii}))
end

cellTypeNormTC = cell(5,1);
cellTypeAvgTC = cell(5,1);
for jj = 1:5
    cellTypeNormTC{jj} = cell(length(cellTypeIndices{jj}),1);
    cellTypeAvgTC{jj} = zeros(size(datarun.stas.time_courses{1}));
    %onMidgNorm = cell(length(iOnM),1); %normalized time courses
    %onMidgAvg = zeros(size(datarun.stas.time_courses{1}));
    for ii = 1:length(cellTypeIndices{jj})
        cellTypeNormTC{jj}{ii} = zeros(size(datarun.stas.time_courses{cellTypeIndices{jj}(ii)}));
        tcSum = sum(sum(abs(datarun.stas.time_courses{cellTypeIndices{jj}(ii)})));
        cellTypeNormTC{jj}{ii} = datarun.stas.time_courses{cellTypeIndices{jj}(ii)}/tcSum;
        cellTypeAvgTC{jj} = cellTypeAvgTC{jj} + cellTypeNormTC{jj}{ii};
    end
    cellTypeAvgTC{jj} = cellTypeAvgTC{jj}/length(cellTypeIndices{jj});
    
    
    figure('position', [100 100 160*1.618 200])
    axes('units', 'points', 'position', [30 30 100*1.618 100])
    hold on
    plot(-8.3*23:8.3:0, cellTypeAvgTC{jj}(:,1), 'r', 'LineWidth', 2)
    plot(-8.3*23:8.3:0, cellTypeAvgTC{jj}(:,2), 'g', 'LineWidth', 2)
    plot(-8.3*23:8.3:0, cellTypeAvgTC{jj}(:,3), 'b', 'LineWidth', 2)
    hold off
    title([datarun.cell_types{jj}.name 'time course'])
    xlabel('ms')
end


% calculates average rf diameter for each cell type

%factor of 10 in area calculation because stixel size = 10 pixels

%gets rid of crappy SBC
for ii = length(cellTypeIndices{5}):-1:1
    if datarun.cell_ids(cellTypeIndices{5}(ii)) == 995
        cellTypeIndices{5}(ii) = [];
    end
end


avgDiamCellTypes = zeros(5,1);
for jj = 1:5
    for ii = 1:length(cellTypeIndices{jj})
        avgDiamCellTypes(jj) = avgDiamCellTypes(jj) + 2*10*sqrt(prod(datarun.stas.fits{cellTypeIndices{jj}(ii)}.sd));
    end
    avgDiamCellTypes(jj) = sqrt(areaRatio)*avgDiamCellTypes(jj)/length(cellTypeIndices{jj});
end

%% load dataset used for center-of-mass calculation


datarun = load_data('2008-08-27-5/data003/data003/data003');
%datarun = load_data('2008-08-27-5', 'rf-4-gf');
datarun = load_index(datarun);

datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun, 'load_sta', [], 'save_rf', 1);
datarun = load_ei(datarun, []);
datarun = load_neurons(datarun);

datarun = get_sta_fits_from_vision(datarun,'all');


cellTypeInd = cell(5,1);
for ii = 1:5
    cellTypeInd{ii} = get_cell_indices(datarun, {ii});
end

%% region defining sub-set of electrodes

arrayPoly = [195 315 195 -45 -165 -45;
            -330 -90 150 150 -90 -330];
        
% shifts region slightly to include an additional SBC, but edge might not have as complete of mosaics
arrayPoly(1,:) = arrayPoly(1,:) + 30;
arrayPoly(2,:) = arrayPoly(2,:) + 60;


        
%% soma location figure

comAvgs = cell(5, 1);

lineColors(1,:) = [90 156 0]/255; %pale grass
lineColors(2,:) = [255 124 59]/255; %salmon
lineColors(3,:) = [101 52 255]/255; %purple
lineColors(4,:) = [52 198 247]/255; %aqua
lineColors(5,:) = [238 55 128]/255; %calm magenta

figure
hold on
plot(datarun.ei.position(:,1), datarun.ei.position(:,2), 'k.')
plot([arrayPoly(1,:) arrayPoly(1,1)], [arrayPoly(2,:) arrayPoly(2,1)], 'k-')

cellsInRegion = cell(5, 1);

for jj = 1:5
    cellIDs = datarun.cell_ids(cellTypeInd{jj});
    %cellIDs = datarun.cell_types{jj}.cell_ids;
    cellType = datarun.cell_types{jj}.name;
    
    cellsInRegion{jj} = false(length(cellIDs), 1);


    comAll = cell(length(cellIDs), 1);
    dists = zeros(length(cellIDs), 1);
    comAvgs{jj} = zeros(2, length(cellIDs));
    for ii = 1:length(cellIDs)
        %figure('position', [200 200 600 600])
        com = ei_com_(get_ei(datarun, cellIDs(ii)), datarun.ei.position, datarun.ei.nlPoints, 'roi', {'peak', 90}, 'frames', [0 20]);
        comAll{ii}(1,:) = com{1}(1,:); %COM of negative portion at frame 0
        comAll{ii}(2,:) = com{2}(2,:); %COM of positive portion at frame 20

        dists(ii) = norm(comAll{ii}(1,:) - comAll{ii}(2,:));
        comAvgs{jj}(1, ii) = mean([comAll{ii}(1,1) comAll{ii}(2,1)]);
        comAvgs{jj}(2, ii) = mean([comAll{ii}(1,2) comAll{ii}(2,2)]);
        
        if any(isnan(comAvgs{jj}(:,ii)))
            warning(['no com calculated for cell' num2str(cellIDs(ii))]) %#ok<WNTAG>
        end
    end

    firstCell = 1;
    for ii = 1:length(cellIDs)
        if inpolygon(comAvgs{jj}(1, ii), comAvgs{jj}(2, ii), arrayPoly(1,:), arrayPoly(2,:))
            cellsInRegion{jj}(ii) = true;
            plot(comAvgs{jj}(1, ii), comAvgs{jj}(2, ii), '.', 'MarkerFaceColor', lineColors(jj,:), 'MarkerEdgeColor', lineColors(jj,:))
            
            %plots circle showing average RF diameter
%             if firstCell
%                 nPoints = 100;
%                 theta = linspace(0, 2*pi, nPoints);
%                 r = ones(1,nPoints)*avgDiamCellTypes(jj)/2;
%                 [x,y] = pol2cart(theta, r);
%                 x = x + comAvgs{jj}(1, ii);
%                 y = y + comAvgs{jj}(2, ii);
%                 plot(x,y, 'color', lineColors(jj,:))
%                 firstCell = 0;
%             end
        end
    end
end

hold off
axis equal

%% lasso particular cell locations to get cell ids and then look at eis
ii = 3;

pointsToPlot = comAvgs{ii}(:, cellsInRegion{ii});
cellIDs = datarun.cell_ids(cellTypeInd{ii});
cellIDs = cellIDs(cellsInRegion{ii}); %removes cells not within array boundary
figure
hold on
plot(datarun.ei.position(:,1), datarun.ei.position(:,2), 'm.')
iPoints = lassoPoints(pointsToPlot'); %returns indices of selected cells
hold off
close

disp(num2str(cellIDs(iPoints)))

plot_ei_scroll(datarun, 604)
plot_ei_scroll(datarun, 873)
com1_0 = ei_com_(get_ei(datarun, 604), datarun.ei.position, datarun.ei.nlPoints, 'roi', {'peak', 90}, 'frames', 0, 'foa', 0);
com1_20 = ei_com_(get_ei(datarun, 604), datarun.ei.position, datarun.ei.nlPoints, 'roi', {'peak', 90}, 'frames', 20, 'foa', 0);
com2_0 = ei_com_(get_ei(datarun, 873), datarun.ei.position, datarun.ei.nlPoints, 'roi', {'peak', 90}, 'frames', 0, 'foa', 0);
com2_20 = ei_com_(get_ei(datarun, 873), datarun.ei.position, datarun.ei.nlPoints, 'roi', {'peak', 90}, 'frames', 20, 'foa', 0);

com1_0 = com1_0{1}; %just keep negative COM
com1_20 = com1_20{2}; %just keep positive COM
com2_0 = com2_0{1};
com2_20 = com2_20{2};

figure
hold on
plot(datarun.ei.position(:,1), datarun.ei.position(:,2), 'k.')
plot([com1_0(1) com1_20(1)], [com1_0(2) com1_20(2)], 'r-')
plot([com2_0(1) com2_20(1)], [com2_0(2) com2_20(2)], 'r-')
hold off

%% plot mosaics of included cells

figure
hold on
for jj = 1:5
    figure
    plot_rf_summaries(datarun, datarun.cell_ids(cellTypeInd{jj}(cellsInRegion{jj})), 'fit_color', lineColors(jj,:), 'plot_fits', true)
    title('data004')
end
hold off


% figure
% hold on
% for jj = 1:5
%     plot_rf_summaries(datarun, datarun.cell_ids(cellTypeInd{jj}), 'fit_color', lineColors(jj,:), 'plot_fits', true)
% end
% hold off


% figure
% hold on
% for jj = 1:5
%     plot(0.1, 0.1*(6-jj), 's', 'MarkerFaceColor', lineColors(jj,:), 'MarkerEdgeColor', lineColors(jj,:))
%     text(0.2, 0.1*(6-jj), datarun.cell_types{jj}.name)
% end
% hold off
% set(gca, 'ylim', [0 1])


%% plots classification scatter plots (unfinished!)


% %code from MARTIN
% paramsFile = edu.ucsc.neurobiology.vision.io.ParametersFile(datarun.names.rrs_params_path);
% 
% for i=1:length(datarun.cell_ids)
%     t = paramsFile.getArrayCell(datarun{j}.cell_ids(i), 'GreenTimeCourse');
%     t1 = paramsFile.getArrayCell(datarun{j}.cell_ids(i), 'RedTimeCourse');
%     t2 = paramsFile.getArrayCell(datarun{j}.cell_ids(i), 'BlueTimeCourse');
%     
%     datarun.green_time_course{i} = t/sqrt(sum(t.^2));
%     t = paramsFile.getArrayCell(datarun{j}.cell_ids(i), 'Auto');
%     datarun.auto{i} = t/sqrt(sum(t.^2));
% 
%     if 1
%         [junk t]=max(max(abs(datarun{j}.ei.eis{i}),[],2));
%         tt=datarun{j}.ei.eis{i}(t,:);
%     else
%         [junk t]=sort(max(abs(datarun{j}.ei.eis{i}),[],2));
%         tt=mean(datarun{j}.ei.eis{i}(t(end-2:end),:));
%     end
% 
%     tt=tt/sqrt(sum(tt.^2));
%     [junk t]=min(tt);
%     datarun{j}.ei_waveform{i}=tt(t-8:t+37)';
% end
% 
%            
% allCellsInRegion = [];
% 
% for ii = 1:5
%     allCellsInRegion = [allCellsInRegion, cellTypeInd{ii}(cellsInRegion{ii})];
% end
% 
% datarun = get_sta_summaries(datarun, 'all');
% 
% plot_rf_summaries(datarun, datarun.cell_ids(allCellsInRegion), 'fit_color', lineColors(jj,:), 'plot_fits', true)
% 
% 
% timeCourses = 
% 
% cellTypeNormTC{jj}{ii} = zeros(size(datarun.stas.time_courses{cellTypeIndices{jj}(ii)}));




