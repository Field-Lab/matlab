% compute and plot contours for a mosaic of RFs

data_spec = '2009-02-28-2/data010/data010';
cell_spec = {3};

if ~exist('datarun','var')  % load data
    
    % initialize struct
    datarun = load_data(data_spec);
    
    % load STAs
    datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
    
    % load cell types
    datarun = load_params(datarun,'verbose',1);
    
    % compute RFs, etc
    datarun = get_sta_summaries(datarun,cell_spec);
    
    % filter summary frames (saved to datarun.stas.summaries_filt{})
    datarun = get_rfs_filtered(datarun,cell_spec,'verbose',1,'filt_params',struct('filt_type','gauss','radius',.4));
    
    % Do Delaunay Triangulation; for determining RoI for calculating UI
    datarun = do_delaunay_tri(datarun, cell_spec);
    
    return;
end

% Pick up where we left off...

% Do Delaunay Triangulation; for determining RoI for calculating UI
datarun = do_delaunay_tri(datarun, cell_spec);

switch 2
    case 1
        % Cut triangles with excessively long edges; indicates where cells are
        % likely missing in mosaic.
        datarun = cull_delaunay_tri(datarun, cell_spec, 1.9);
        
        % Build Region of Interest from culled Delaunay Triangulation
        datarun = build_rf_roi(datarun, cell_spec);
        
        % Get contour threshold that maximizes uniformity index
        datarun = maximize_uniformity_index(datarun, cell_spec, 0.3, 0.05, 0.01);
        
        % compute contours (summary frames normalized, then saved to datarun.stas.contours)
        datarun = get_rf_contours(datarun, cell_spec, datarun.stas.mosaics{cell_spec{:}}.best_thresh, 'norm_params', struct('method', 'peak'), 'verbose', 1);
        
        % plot contours
        plot_rf_summaries(datarun, cell_spec, 'plot_contours', 1, 'foa', 1, 'label', 0);
        lock;
        
        % plot simplified contours, with alpha fills
        datarun = simplify_rf_contours(datarun, cell_spec);
        plot_rf_summaries(datarun, cell_spec, 'plot_contours', 1, 'contours_field', 'rf_contours_simple', 'foa', 1, 'label', 0, 'contour_fill', 'r', 'contour_alpha', 0.5);
    case 2
        % ALTERNATIVELY...
        plot_rf_coloring(datarun, cell_spec, 'rfs', 'summaries_filt');
end
