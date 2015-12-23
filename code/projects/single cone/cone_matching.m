%%

piece = '2012-01-27-1';
dataruns = {'d1' 'd2'};
d1 = load_data(['/snle/acquisition/' piece '/data002/data002']);
d2 = load_data(['/snle/acquisition/' piece '/data004/data004']);

for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    datarun = load_params(datarun);
    datarun = load_neurons(datarun);
    assignin('caller', dataruns{i}, datarun);
end

for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    datarun = load_ei(datarun, 'all', 'keep_java_ei', false);
    assignin('caller', dataruns{i}, datarun);
end

for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});
    
    datarun = load_sta(datarun,'load_sta',[], 'keep_java_sta', false, 'guess_stimulus', false);
    datarun = set_polarities(datarun);
    datarun = get_sta_fits_from_vision(datarun);
    
    assignin('caller', dataruns{i}, datarun);
end


% Normal cone loading
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});

    datarun = load_cones(datarun, 1);
    
    assignin('caller', dataruns{i}, datarun);
end


% Mosaic generation
for i = 1:length(dataruns)
    datarun = evalin('caller', dataruns{i});

    datarun = make_mosaic_struct(datarun);
    datarun = make_voronoi_masks(datarun);
    
    assignin('caller', dataruns{i}, datarun);
end

clear dataruns datarun i piece


%% Cell picks
onM_d00  = [1561 1608 1817 1982 2283 2612 2658 2716 2822 3016 3196 3617 3961 4517 5191 6392];
offM_d00 = [3937 812 1006 707 5419];

onM_d04 = [601 1562 1681 1968 2283 2658 2821 3244 4562 6542 6601 7336];
offM_d04 = [541 5161];


%%
rgc_match = map_ei(d1, d2);


%% Cone match
cell_types = {{4}};
arrow_colors = {'w' 'k'};

overlayf = figure();
overlay_cone_mosaics(d1, d2);
hold on;

cmatches = cell(size(d1.cones.centers,1), 1);
for i = 1:length(cell_types)
    cell_type = cell_types{i};
    arrow_color = arrow_colors{i};

    rgcs = d1.cell_ids(get_cell_indices(d1, cell_type));
    xcshifts = rf_xcorr_shifts(d1, d2, rgcs, 'interp_delta', 1);
    
    for j = 1:length(rgcs)
        rgc = rgcs(j);
        xcshift = xcshifts(j,:);
        if isnan(xcshift(1)), continue; end
        
        [~, selection] = select_cone_weights(d1, rgc, 'thresh', 0.25, 'radius', [0 inf], 'polarity', 1, 'contiguity', true);
%         cones1 = find(selection);
%         neighbors = ipdm(d1.cones.centers(cones1,:)+repmat(xcshift, [length(cones1) 1]), d2.cones.centers, 'Subset', 'NearestNeighbor', 'Result', 'Structure');
        
        figure(overlayf);
%         for k = 1:length(cones1)
%             dist = neighbors.distance(k);
%             if dist > 2, continue; end
%             
%             r1 = neighbors.rowindex(k);
%             cone1 = cones1(r1);
%             cmatch = neighbors.columnindex(k);
%             cmatches{cone1}(end+1) = cmatch;
%             plot([d1.cones.centers(cone1,1) d2.cones.centers(cmatch,1)], [d1.cones.centers(cone1,2) d2.cones.centers(cmatch,2)], arrow_color);
%         end
        quiver(d1.cones.centers(selection,1), d1.cones.centers(selection,2), repmat(xcshift(1), [sum(selection) 1]), repmat(xcshift(2), [sum(selection) 1]), 0, arrow_color);
        drawnow;
        
    end; clear rgc match ax1 ax2
end; clear i j k

clear cell_type cell_types


%% Get cone closest to click
[x,y] = ginput(1);
near1 = ipdm([x y], d1.cones.centers, 'Subset', 'NearestNeighbor', 'Result', 'Structure');
near2 = ipdm([x y], d2.cones.centers, 'Subset', 'NearestNeighbor', 'Result', 'Structure');
if near1.distance < near2.distance
    cn = near1.columnindex;
    title(['Map1 ' num2str(cn)]);
else
    cn = near2.columnindex;
    title(['Map2 ' num2str(cn)]);
end


