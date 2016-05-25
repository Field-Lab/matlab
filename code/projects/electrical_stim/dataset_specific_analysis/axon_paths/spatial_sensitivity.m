% Comparing distance of electrode from predicted axon path and activation
% threshold

% % cells ids from /Volumes/Analysis/2015-04-14-0/data000/data000.ei
% cellIds = [1 392 633 917 1308 1699 2345 3049 3466 3947 4655 6873];
[xc,yc] = getElectrodeCoords512(); 
codebase_path = matlab_code_path; 
cmap = load([codebase_path 'code/projects/electrical_stim/'...
    'resources/redtealcmap.mat']); 
% cells ids from /Volumes/Analysis/2012-09-24-3/data000/data000.ei
dataPath = '/Volumes/Analysis/2012-09-24-3/data000/data000'; 
cellIds = [3457 4058 4656 2796 5123 5748 6439];

% Analyze all electrodes that recorded voltages over some threshold. 
% Load EI. 
datarun  = load_data(dataPath);
datarun  = load_neurons(datarun);
datarun  = load_sta(datarun, 'load_sta', 'all');
datarun  = load_params(datarun);
datarun  = load_ei(datarun, cellIds,'keep_java_ei','false');

x_coords = [];
y_coords = [];
stim_elecs = [];

ei_thresh = 10; 
%%
figure(100); 
%figure(200); scatter(xc,yc,3,'k','filled');    axis off; axis image;
%line([0 0],[min(yc) max(yc)],'Color','k');
%line([min(xc)/2 min(xc)/2],[min(yc) max(yc)],'Color','k');
%line([max(xc)/2 max(xc)/2],[min(yc) max(yc)],'Color','k');

distance_threshold = [];

for n = 1:length(cellIds)
    thresh_quad1 = zeros(1,512);
    thresh_quad2 = zeros(1,512);
    thresh_quad3 = zeros(1,512);
    thresh_quad4 = zeros(1,512);
    neuronId = cellIds(n);
    ei = datarun.ei.eis{get_cell_indices(datarun,neuronId)};
    eiAmps = max(ei,[],2) - min(ei,[],2);
    if n == 1
        plotCoords = true;
    else
        plotCoords = false;
    end; 
    [ax, ay] = eiContour_wPolyFit(eiAmps,'figureNum',100,'ei_thresh', ei_thresh,'plotCoords',plotCoords);
    
    % Find ei amps greater than a particular threshold
    elecs = find(eiAmps > ei_thresh); 
    %figure(200); 
    %hold on; scatter(xc(elecs),yc(elecs),eiAmps(elecs),'r','filled'); 
 
    % Determine the quadrant for each electrode in a given cell .  
    
    if ~isempty(find(xc(elecs) < min(xc/2), 1))
        fpath = '/Volumes/Analysis/2012-09-24-3/data003/';
        [~, thresh_quad1] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',0,'suppressPlots', 1);
    end

    if ~isempty(find(xc(elecs) > min(xc/2) & xc(elecs) < 0, 1));
        fpath = '/Volumes/Analysis/2012-09-24-3/data004/';
        [~, thresh_quad2] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',0,'suppressPlots', 1);
    end
    
    if ~isempty(find(xc(elecs) < max(xc/2) & xc(elecs) > 0, 1))
        fpath = '/Volumes/Analysis/2012-09-24-3/data005/';
        [~, thresh_quad3] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',0, 'suppressPlots', 1);

    end
    if ~isempty(find(xc(elecs) > max(xc/2), 1))
        fpath = '/Volumes/Analysis/2012-09-24-3/data006/';
        [~, thresh_quad4] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',0, 'suppressPlots', 1);

    end
   
    thresholds = thresh_quad1 + thresh_quad2 + thresh_quad3 + thresh_quad4;
    
    thresholds(thresholds>6) = 0; 
    thresholds(thresholds<1) = 0; 
    idx = find(thresholds);
    
    
    for elec = 1:length(idx)
        min_dist_to_axon = sqrt( min( sum( bsxfun(@minus, vertcat(ax,ay), [xc(idx(elec));yc(idx(elec))]).^2,1)));
        point = [min_dist_to_axon thresholds(idx(elec))];
        distance_threshold = vertcat(distance_threshold, point);
    end
        
    figure(100); 
    
    x_coords = [x_coords xc(idx)];
    y_coords = [y_coords yc(idx)];
    stim_elecs = [stim_elecs thresholds(idx)];
    
    hold on; 
  
    
end

scatter(x_coords,y_coords,200, stim_elecs,'filled');
colormap(cmap.m); caxis([0 5]);
colorbar;



figure(300);
scatter(distance_threshold(:,1), distance_threshold(:,2));
xlabel('Distance between stimulus electrode and nearest point on predicted axon path (um)');
ylabel('Axon stimulation threshold (uA)');

