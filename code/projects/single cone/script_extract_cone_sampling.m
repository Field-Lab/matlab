% script to identify the locations and spectral properties of cones in fine RFs
%
%
% procedure
%
% * look across all cells to get cone locations
% * at each cone, figure out if it's L, M, or S
% * explain each RF as a sum of individual cones
% * make nice plots and save data in text files
%
%
%  all parameters are at the top of the script
%  options to execute each component are scattered throughout the script
%
%
% gauthier 2008-10
%



% continue when there is an error?
%catch_errors = true; % the error is printed and the code continues to the next datarun
catch_errors = false; % everything stops

% run 1003 later, had bug when plotting cell sampling
 dataset_list = 1007;

for dd = dataset_list
    
    % save_suffix
    clear save_suffix

    save_suffix = '-old_corrected';

    try


%%     	SET PARAMETERS



        %%%%%%%%%%%%%%%   loading/saving datarun   %%%%%%%%%%%%%%%

        % where datarun .mat files are stored
        %save_datarun_prefix = '/snle/home/gauthier2/Desktop/2008-one/saved/save-';
        save_datarun_prefix = '/Users/gfield/Desktop';
        %save_datarun_prefix = [single_cone_path 'save-']; % commented out
        %by gdf

        % load data from vision files rather than looking for a mat file?
        load_fresh = true; % load from vision files "from scratch"
        %load_fresh = false; % look for a mat file (if none is found, load from vision files)


        % these arguments only apply if load_fresh is true or no mat file is found:

        % parameters to load marks
        clear marks_params
        marks_params.strength = {'inner or',[0.4044 0.8854 0.2292;0.1483 0.9161 0.3726;0.0174 0.0792 0.9967]};
            
        % keep all RFs?  if false, all other summary info (e.g. significant stixels) will still be loaded
        keep_rfs = true;
        %keep_rfs = false;
        
        % save datarun to disk?
        save_datarun = true;
        %save_datarun = false;
        
        % print portraits?
        save_portraits = true;
        %save_portraits = false;

        % load special set of cells
        % parasol, midget, and SBC by default.  can be changed for a specific dataset below.
        %load_cell_spec = {1,2,3,4,5};
        load_cell_spec = 'all';




        %%%%%%%%%%%%%%%   saving results   %%%%%%%%%%%%%%%

        % where to save files (directory with trailing slash)
        save_file_dir = single_cone_path;
        
        % write text files?
        %save_results = true;
        save_results = false;

        % prompt the user to select ROI for RGCs
        export_rgc_roi = false;
        %export_rgc_roi = true;





        %%%%%%%%%%%%%%%   which cells to analyze   %%%%%%%%%%%%%%%

        % which cell types to use to get cone locations and colors
        %cone_finding_cell_spec = {1,2,3,4,5};
        cone_finding_cell_spec = 'all';

        % which cell types to use to identify the sampling of
        %regression_cell_spec = {1,2,3,4,5};
        regression_cell_spec = 'all';


        % regions of interest for cells and cones
        %
        % set a value to 'roi' to use only a region of interest.  the user will be prompted
        % to click points on an image to define a ROI.

        which_cells = 'all';
        %which_cells = 'roi';

        which_cones = 'all';
        %which_cones = 'roi';


        
        %%%%%%%%%%%%%%%   blur RFs   %%%%%%%%%%%%%%%
        
        %blur_rfs = true;
        blur_rfs = false;
        
        % how much to blur by
        blur_radius = 1.3;
        
        % plot before and after
        blur_figure = 89;




        %%%%%%%%%%%%%%%   spatial sensitivity combined across all RGCs   %%%%%%%%%%%%%%%

        clear spat_sens_params

        % how to combine RGB values of the RF
        spat_sens_params.strength = {'inner or',...
            [0.4044 0.8854 0.2292;0.1483 0.9161 0.3726;0.0174 0.0792 0.9967]};

        % how to filter the RF before looking for significant stixels
        spat_sens_params.filter = [];

        % how to find significant stixels
        spat_sens_params.selection_params = struct('type','thresh','thresh',5);
        %params.selection_params = struct('type','max','thresh',5);

        % how to combine stixels between cells
        spat_sens_params.combine_stixels = 'sum';
        %cone_loc_params.combine_stixels = 'max';

        % online readout of what's going on
        spat_sens_params.verbose = true;
        spat_sens_params.fig_single_cell = [];
        spat_sens_params.foa_spat_sens = 102;


        
        %%%%%%%%%%%%%%%   how to combine cones across RFs   %%%%%%%%%%%%%%%

        % how to combine cones across RFs
        combine_rfs = 'cluster'; % find cones in each RF, then cluster locations
        %combine_rfs = 'spat sens'; % compute spatial senstivity
        
        % this params struct is only used if combine_rfs == 'cluster'
        clear cluster_rf_params
        cluster_rf_params.radius_1 = 1.2;
        cluster_rf_params.radius_2 = 1.6;
        


        %%%%%%%%%%%%%%%   cone number and location seeds   %%%%%%%%%%%%%%%

        % how to find cone center points
        cone_finding = 'maxima';
        %cone_finding = 'segment';




        %%%%%%%%%%%%%%%   cone location and color identification   %%%%%%%%%%%%%%%

        clear cone_rf_params cone_kernel_params

        % how to identify cone center points
        %cone_rf_params.centers = 'fit'; % fit a gaussian
        %cone_rf_params.centers = 'com'; % find the center of mass
        %cone_rf_params.centers = 'fixed';  % use the initial conditions
        cone_rf_params.centers = 'fit neighborhood';  % fit all cones in a small neighborhood
        
        % how to combine RGB values of the RF
        cone_rf_params.sensitivity.strength = {'inner or',...
            [0.4044 0.8854 0.2292;0.1483 0.9161 0.3726;0.0174 0.0792 0.9967]};

        % how to filter the RF before looking for significant stixels
        cone_rf_params.sensitivity.filter = [];
        

        % cone kernel parameters
        % this shape is fitted to each cone
        cone_kernel_params.type = 'dog';
        cone_kernel_params.center_radius = 0.75;
        cone_kernel_params.surround_scale = 0;
        cone_kernel_params.surround_radius = 1;

        % how to extract a single RGB triplet from several RGCs
        % (see summarize_cone_rfs and rbg_from_cones)
        switch 2
            case 1 % PCA on all stixels
                cone_rf_params.regress = 0;
                cone_rf_params.combine = 'pca';

            case 2 % sum based on regression values
                cone_rf_params.regress = 1;
                cone_rf_params.combine = 'sum';
        end

        % online readout of what's going on
        cone_rf_params.single_cone_figure = [];%105;
        cone_rf_params.cone_weights_figure = [];%106;
        cone_rf_params.verbose = true;




        %%%%%%%%%%%%%%%   classify cones   %%%%%%%%%%%%%%%

        clear classify_params

        classify_params.algorithm = 'k means, EM';




        %%%%%%%%%%%%%%%   plot cone mosaic   %%%%%%%%%%%%%%%

        clear plot_mosaic_params

        plot_mosaic_params.fig_or_axes = 108;





        %%%%%%%%%%%%%%%   visualize cone projections   %%%%%%%%%%%%%%%

        % which 2D projection space to use to visualize cone colors

        clear remap
        switch 2
            case 1 % R/(R+G)  vs  B/(B+G)
                remap.fcn = @(x)([x(:,1)./(x(:,1)+x(:,2)) x(:,3)./(x(:,3)+x(:,2))]);
                remap.x_caption = 'red/(red+green)';
                remap.y_caption = 'blue/(blue+green)';
            case 2 % R/G  vs  B/G
                remap.fcn = @(x)([x(:,1)./x(:,2) x(:,3)./x(:,2)]);
                remap.x_caption = 'red/green';
                remap.y_caption = 'blue/green';
        end











%%      LOAD DATARUN




        % clear previous datarun, closing the java_sta object first
        if exist('datarun','var')
            if isfield(datarun,'stas') && isfield(datarun.stas,'java_sta') && isjava(datarun.stas.java_sta)
                datarun.stas.java_sta.close;
            end
            clear datarun
        end
            
        % clear previous short name
        clear short_name

        clear rig light_path
        switch dd
            
            
            % DATASETS FOR PAPER
            
            
            case 1001 % plantain
                dataset_mat_file = [single_cone_path 'saved/plantain.mat'];
                dataset_spec = {'plantain'};
            
            case 1002 % peach
                dataset_mat_file = [single_cone_path 'saved/peach.mat'];
                dataset_spec = {'peach'};
            
            case 1003 % blueberry
                dataset_mat_file = [single_cone_path 'saved/blueberry.mat'];
                dataset_spec = {'blueberry'};
            
            case 1004 % apricot
                dataset_mat_file = [single_cone_path 'saved/apricot.mat'];
                dataset_spec = {'apricot'};
            
            case 1005 % kiwi
                dataset_mat_file = [single_cone_path 'saved/kiwi.mat'];
                dataset_spec = {'kiwi'};
            
            case 1006 % grapes
                dataset_mat_file = [single_cone_path 'saved/grapes.mat'];
                dataset_spec = {'grapes'};
            
            case 1007 % apple
                dataset_mat_file = [single_cone_path 'saved/apple-13.mat'];
                dataset_spec = {'apple'};
            
            
                
                
            
            
            case 101 % apricot
                dataset_spec = {'2009-04-13-5/data005-3600s-7200s/data005/data005'};
                rig = 'A'; light_path = 'below';
                
            case 102 % blueberry
                dataset_spec = {'2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551'};
                rig = 'A'; light_path = 'below';
                
                error('already analyzed')
                
            case 103 % butterfly
                dataset_spec = {'2008-12-12-1/data005/data005-3600s-7200s/data005/data005'};
                rig = 'A'; light_path = 'below';
                
            case 104 % cherimoya
                dataset_spec = {'2008-03-25-3/data002/data002-0s-3000s/data002/data002'};
                rig = 'A'; light_path = 'above';
                
            case 105 % cherry
                dataset_spec = {'2009-02-28-0/data006/data006-0s-3600s/data006/data006'};
                rig = 'B'; light_path = 'below';
                
            case 106 % grapes
                dataset_spec = {'2007-03-27-2/data014/data014-6180s-13380s/data014/data014'};
                rig = 'B'; light_path = 'below';
                
            case 107 % kiwi
                dataset_spec = {'2008-05-13-3/data006-3600s-7200s/data006/data006'};
                rig = 'B'; light_path = 'below';
                
            case 108 % mango
                dataset_spec = {'2008-04-30-2/data004/data004-0s-3600s/data004/data004'};
                rig = 'A'; light_path = 'below';
                spat_sens_params.filter = struct('type','given','filt',...
                    make_gaussian('center_radius',0.6,'x_size',5,'y_size',5,'center',[3 3]));
                
            case 109 % peach, last 60
                %dataset_spec = {'2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260'};
                rig = 'A'; light_path = 'below';
                
            case 110 % peach, first 40
                dataset_spec = {'2008-08-27-0/data001-0s-2400s/data001/data001'};
                rig = 'A'; light_path = 'below';
                
            case 111 % plantain
                dataset_spec = {'2008-08-27-5/data003/data003-3600s-7200s/data003/data003'};
                rig = 'A'; light_path = 'below';
                
            case 112 % plum
                dataset_spec = {'2008-04-22-5/data006-1800s-3600s/data006/data006'};
                rig = 'B'; light_path = 'below';
                
            case 113 % pomegranate
                %dataset_spec = {'2007-08-21-1/data003/data003'};
                rig = 'B'; light_path = 'below';
            
            
                
                
                
                
                
            
            
            case 1 % plantain
                dataset_spec = {'2008-08-27-5/data003/data003/data003'};
                rig = 'A'; light_path = 'below';
                load_cell_spec = 'all';
                cone_finding_cell_spec = 'all';
                regression_cell_spec = 'all';
            case 2 % grapes (from faseb 2008 poster)
                dataset_spec = {'2007-03-27-2/data014/data014/data014'};
                rig = 'B'; light_path = 'below';
                %spat_sens_params.filter = struct('type','given','filt',...
                %    make_gaussian('center_radius',0.6,'x_size',5,'y_size',5,'center',[3 3]));
            case 3 % kiwi
                dataset_spec = {'2008-05-13-3/data006/data006'};
                rig = 'B'; light_path = 'below';
            case 4 % mango
                dataset_spec = {'2008-04-30-2/data004/data004/data004'};
                rig = 'A'; light_path = 'below';
                spat_sens_params.filter = struct('type','given','filt',...
                    make_gaussian('center_radius',0.6,'x_size',5,'y_size',5,'center',[3 3]));
                load_cell_spec = {1,2,3,4,5};
                cone_finding_cell_spec = 'all';
                regression_cell_spec = 'all';
                % to fix cell types:
                %onp = [1356 738 1263 407 7084 6722 7043 6273 5254]
                %datarun.cell_types{1}.cell_ids = sort([datarun.cell_types{1}.cell_ids onp]);
                %datarun.cell_types{3}.cell_ids = setdiff(datarun.cell_types{3}.cell_ids,onp);

            case 5 % pomegranate
                dataset_spec = {'2007-08-21-1/data003/data003'};
                rig = 'B'; light_path = 'below';
            case 6 % plum
                dataset_spec = {'2008-04-22-5/data006/data006'};% mem
                rig = 'B'; light_path = 'below';
                %load_cell_spec = {2,4,5,8,9,10};
                %cone_finding_cell_spec = load_cell_spec;
                %regression_cell_spec = load_cell_spec;
            case 7 % cherimoya
                dataset_spec = {'2008-03-25-3/data002/data002'};
                rig = 'A'; light_path = 'above';
            case 8 % raspberry, >6000 cones
                %dataset_spec = {'2006-06-12-0','rf-3-auto'}; % pca
                dataset_spec = {'2006-06-12-0/data003/data003'}; % nwpca
                rig = 'B'; light_path = 'above';
                
                
            case 27 % kiwi with pca
                dataset_spec = {'2008-05-13-3/data006-pca/data006'};
                rig = 'B'; light_path = 'below';
                
            case 28 % peach second hour
                dataset_spec = {'2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260'};
                rig = 'A'; light_path = 'below';
                
            case 29 % butterfly
                dataset_spec = {'2008-12-12-1/data005/data005'};
                rig = 'A'; light_path = 'below';
                
            case 30 % 2008-12-12-0
                dataset_spec = {'2008-12-12-0/data004/data004'};
                rig = 'B'; light_path = 'above';
                
            case 31 % cherry
                dataset_spec = {'2009-02-28-0/data006/data006'};
                rig = 'B'; light_path = 'below';
                load_cell_spec = {1,2,3,4,5};
                cone_finding_cell_spec = {1,2,3,4,5};
                regression_cell_spec = {1,2,3,4,5};
                
                
            case 47 % unnamed
                dataset_spec = {'2009-04-13-1/data004/data004'};
                rig = 'B'; light_path = 'below';
                
            case 48 % unnamed, 1st hour
                dataset_spec = {'2009-04-13-1/data004-s0-s3599/data004-s0-s3599'};
                rig = 'B'; light_path = 'below';
                
            case 49 % unnamed, 2nd hour 
                dataset_spec = {'2009-04-13-1/data004-s3600-/data004-s3600-'};
                rig = 'B'; light_path = 'below';
                
            case 50 %
                dataset_spec = {'2009-04-13-4/data005/data005'};
                rig = 'B'; light_path = 'below';
                
            case 51 % apricot
                dataset_spec = {'2009-04-13-5/data005/data005'};
                rig = 'A'; light_path = 'below';
                
            case 52 % apricot, 2nd hour
                dataset_spec = {'2009-04-13-5/data005-s3600-/data005-s3600-'};
                rig = 'A'; light_path = 'below';
                
                
                
            case 32 % peach, second hour
                dataset_spec = {'2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260'};
                rig = 'A'; light_path = 'below';
                cone_finding = 'manual';
                cone_marks = getfield(load([single_cone_path...
                    '2008-08-27-0_data001-s5661-s9260_data001-s5661-s9260-manual1/cone_manual_peach_1s.mat']),...
                    'cone_marks');
                cone_finding_cell_spec = {1,2,3,4,5,6};
                cone_rf_params.centers = 'fixed';
                save_suffix = '-manual1';
                
            case 33 % peach, second hour
                dataset_spec = {'2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260'};
                rig = 'A'; light_path = 'below';
                cone_finding = 'manual';
                cone_marks = getfield(load([single_cone_path...
                    '2008-08-27-0_data001-s5661-s9260_data001-s5661-s9260-manual2/cone_manual_peach_2s.mat']),...
                    'cone_marks');
                cone_finding_cell_spec = {1,2,3,4,5,6};
                load_cell_spec = {1,2,3,4,5,6};
                cone_rf_params.centers = 'fixed';
                save_suffix = '-manual2';
                
                
            case 34 % pomegranate
                dataset_spec = {'2007-08-21-1/data003/data003'};
                rig = 'B'; light_path = 'below';
                cone_finding = 'manual';
                cone_marks = getfield(load([single_cone_path...
                    '2007-08-21-1_data003_data003-manual1/cone_manual_pomergrante_1s.mat']),...
                    'cone_marks');
                cone_finding_cell_spec = {1,2,3,4,5,6};
                cone_rf_params.centers = 'fixed';
                save_suffix = '-manual1';
                
                
            case 35 % pomegranate
                dataset_spec = {'2007-08-21-1/data003/data003'};
                rig = 'B'; light_path = 'below';
                cone_finding = 'manual';
                cone_marks = getfield(load([single_cone_path...
                    '2007-08-21-1_data003_data003-manual2/cone_manual_pomergrante_2s.mat']),...
                    'cone_marks');
                cone_finding_cell_spec = {1,2,3,4,5,6};
                cone_rf_params.centers = 'fixed';
                save_suffix = '-manual2';
                
                
                
            case 36 % plantain, manual "clean"
                dataset_spec = {'2008-08-27-5/data003/data003/data003'};
                rig = 'A'; light_path = 'below';
                cone_finding = 'manual';
                cone_finding_cell_spec = {1,2,3,4,5,6};
                cone_rf_params.centers = 'fixed';
                
                save_suffix = '-manual-clean';
                cone_marks = getfield(load([single_cone_path...
                    '2008-08-27-5_data003_data003_data003' save_suffix '/cone_manual_plantain_clean.mat']),...
                    'cone_marks');
                
                
            case 38 % plantain
                dataset_spec = {'2008-08-27-5/data003/data003/data003'};
                rig = 'A'; light_path = 'below';
                cone_finding = 'manual';
                cone_finding_cell_spec = {1,2,3,4,5,6};
                cone_rf_params.centers = 'fixed';
                
                save_suffix = '-manual-blue';
                cone_marks = getfield(load([single_cone_path...
                    '2008-08-27-5_data003_data003_data003' save_suffix '/cone_manual.mat']),...
                    'cone_marks');
     
                
            case 39 % 
                dataset_spec = {'2008-04-30-1/data004/data004-nwpca/data004'};
                rig = 'B'; light_path = 'below';
                
            case 40 % 
                dataset_spec = {'2008-04-30-1/data004/data004-nwpca/data004'};
                rig = 'B'; light_path = 'below';
                cone_finding = 'manual';
                cone_finding_cell_spec = {1,2,3,4,5,6};
                cone_rf_params.centers = 'fixed';
                
                save_suffix = '-manual-blue';
                cone_marks = getfield(load([single_cone_path...
                    '2008-04-30-1_data004_data004-nwpca_data004' save_suffix '/cone_manual.mat']),...
                    'cone_marks');

                
                
            case 9 % reasonable
                dataset_spec = {'2008-04-22-5/data001/data001'};
                rig = 'B'; light_path = 'below';
            case 10 % ???
                dataset_spec = {'2008-03-25-4/data003/data003'};
                rig = 'B'; light_path = 'below';
            case 11 % ???
                dataset_spec = {'2008-04-30-1/data004/data004'};
                rig = 'B'; light_path = 'below';
            case 12 % ???
                dataset_spec = {'2008-05-13-2/data003/data003'};
                rig = 'A'; light_path = 'above';
            case 13 % ???
                dataset_spec = {'2008-06-10-1/data003-gdf/data003'};
                rig = 'B'; light_path = 'below';
            case 14 % ???
                dataset_spec = {'2008-07-07-3/data005/data005'};
                rig = 'A'; light_path = 'below';

            case 15 % ???
                dataset_spec = {'2006-07-12-0/data010/data010'};
                rig = 'B'; light_path = 'above';
                load_fresh = true;
                save_datarun = false;
                keep_rfs = false;
                
            case 16 % ???
                dataset_spec = {'2006-11-08-1/data003/data003'};
                rig = 'B'; light_path = 'below';
                load_cell_spec = 'all';
                cone_finding_cell_spec = 'all';
                regression_cell_spec = 'all';
            case 17 % ???
                dataset_spec = {'2006-06-15-2/data004/data004'};
                rig = 'B'; light_path = 'above';
            case 18 % ???
                dataset_spec = {'2007-05-01-3/data005/data005'}; % no cell types?
                rig = 'B'; light_path = 'above';
                load_cell_spec = 'all';
                cone_finding_cell_spec = 'all';
                regression_cell_spec = 'all';
            case 19 % ???
                dataset_spec = {'2007-05-01-3/data006/data006'};% not classified
                rig = 'B'; light_path = 'below';
            case 20 % ???
                dataset_spec = {'2007-08-24-1/data010/data010'};
                rig = 'B'; light_path = 'below';

            case 24 % very blurred cones
                dataset_spec = {'2008-08-27-0/data001/data001'};
                rig = 'A'; light_path = 'below';
                
            case 25 % blueberry
                dataset_spec = {'2008-08-26-2/data001/data001'};
                rig = 'A'; light_path = 'below';
                
            case 26 % BW
                dataset_spec = {'2008-08-27-0/data003/data003'};
                rig = 'A'; light_path = 'below';
                
                
                
            case 41 % BW
                dataset_spec = {'2006-07-05-0/data009/data009'};
                rig = 'B'; light_path = 'above';
            case 42 % BW
                dataset_spec = {'2008-04-08-1/data008/data008'};
                rig = 'B'; light_path = 'below';
            case 43 % BW
                dataset_spec = {'2008-05-13-2/data004/data004'};
                rig = 'A'; light_path = 'below';
            case 44 % BW
                dataset_spec = {'2008-05-13-3/data007/data007'};
                rig = 'B'; light_path = 'below';
            case 45 % BW
                dataset_spec = {'2008-06-10-5/data001/data001'};
                rig = 'B'; light_path = 'above';
                
                

                % not on server yet

        end
        
        

        % initialize struct
        datarun = load_data(dataset_spec{:});

        % put short name in more accessible variable
        if ~exist('short_name','var')
            short_name = datarun.names.short_name;
        end

        % display name
        fprintf('\nLoading %s\n',short_name)

        
        
        % identify name of file to load
        if exist('dataset_mat_file','var')
            file_to_load = dataset_mat_file;
        else
            file_to_load = [save_datarun_prefix short_name '.mat'];
        end
        
        
        % if the file exists, load it from disk
        if exist(file_to_load,'file') && ~load_fresh
            
            
            fprintf('Loading %s...\n',file_to_load)

            % load mat file
            load(file_to_load)
            
            % load index file
            datarun = load_index(datarun);
            
            % get java_sta object
            datarun.stas.java_sta = load_java_sta(datarun);

            % show what cell types are present
            show_cell_types(datarun.cell_types)

        else
            % otherwise, get all the data from scratch

            % enter piece and light path info
            if exist('rig','var'); datarun.piece.rig = rig;end
            if exist('light_path','var'); datarun.piece.optical_path_direction = light_path;end

            % fill in stixel height, width
            datarun.stimulus.stixel_height = 1;
            datarun.stimulus.stixel_width = 1;

            % load params file
            datarun = load_params(datarun,struct('verbose',1));
            
            % make class for unclassified cells
            datarun = make_unclassified_cell_type(datarun);

            % load sta file
            datarun = load_sta(datarun,'load_sta',[],'verbose',1,'keep_java_sta',1);
            
            % handle BW STAs
            if strcmpi(datarun.stimulus.independent,'nil')
                marks_params.strength = 'vector length';
            end
            
            
            % load the various STA summaries, don't save STAs
            datarun = get_sta_summaries(datarun,load_cell_spec,'verbose',1,'keep_stas',0,...
                'keep_rfs',keep_rfs,'fig_or_axes',1,'marks_params',marks_params);
            
            % set polarity flags
            datarun = set_polarities(datarun);

            % save to disk
            if save_datarun
                save([save_datarun_prefix short_name],'datarun','-v7.3')
            end
            
            % make summary portraits
            if save_portraits
                for cc =1:6
                    fig_portraits = 2;
                    plot_rf_portraits(datarun,{cc},'figure',fig_portraits,'plot_radius',20,'scale_factor',4)
                    print(fig_portraits,sprintf('%sportraits-%d',save_file_path,cc),'-dpdf')
                end
            end

        end



        % make directory to save in
        
        % if user did not specify a save_suffix, make it empty
        if ~exist('save_suffix','var')
            save_suffix = [];
        end
        
        datarun.names.short_name = short_name_from_names(datarun.names);
        
        % set file path, and create it
        save_file_path = [save_file_dir datarun.names.short_name save_suffix '/'];
        mkdir(save_file_path)





%%   	PERFORM ANALYSIS





        % PLOT EXPECTED CONE SENSITIVITY TO R, G, & B GUNS

        if 0

            % plot expected cone weights
            cone_rgb_expected(datarun,struct('normalize',true,'figure',100));

            % compute remaped weights
            %remap_L = remap.fcn(rgb_weights.L);
            %remap_M = remap.fcn(rgb_weights.M);
            %remap_S = remap.fcn(rgb_weights.S);

        end



        % handle black and white 
        if strcmpi(datarun.stimulus.independent,'nil')

            % expand RFs to occupy all color channels
            for cc =1:length(datarun.cell_ids)
                if ~isempty(datarun.stas.rfs{cc})
                    temp = datarun.stas.rfs{cc};
                    datarun.stas.rfs{cc}(:,:,2) = temp;
                    datarun.stas.rfs{cc}(:,:,3) = temp;
                end
            end
            clear temp
            
            % chose a phony-baloney cone classification algorithm
            classify_params.algorithm='nearest line';
        end
        
        
        
        % BLUR EACH RF
        if blur_rfs
            % define blur kernel
            blur_kernel = fspecial('gaussian',ceil(blur_radius)*3,blur_radius);
            % get cell indices
            cell_indices = get_cell_indices(datarun,'all');
            % for each cell
            for cc = 1:length(cell_indices)
                cell_index = cell_indices(cc);
                % ignore if no RF is present
                if isempty(datarun.stas.rfs{cell_index});continue;end
                
                figure(blur_figure);clf;
                subplot(121);imagesc(norm_image(datarun.stas.rfs{cell_index}));axis image;
                % blur the RF
                datarun.stas.rfs{cell_index} = imfilter(datarun.stas.rfs{cell_index},blur_kernel);
                subplot(122);imagesc(norm_image(datarun.stas.rfs{cell_index}));axis image;
                drawnow
            end
        end
        




        % COMPUTE SPATIAL SENSITIVITY
        % accumulate a set of significant stixels across all cells
        % this is useful to do even if the accumulated spatial sensitivity is not used for cone finding

        switch which_cells
            case 'all'
                % use all cells
                rgc_analysis_roi = [];
            case 'roi'

                % prompt user to select a ROI if it doesn't already exist
                if ~exist('rgc_analysis_roi','var')
                    figure;clf;rgc_analysis_roi = roipoly(zeros(datarun.stimulus.field_height,datarun.stimulus.field_width));
                end
        end

        % use roi
        spat_sens_params.roi = rgc_analysis_roi;

        % identify cone locations for all cells
        [spatial_sensitivity,all_sig_stixels,spatial_cell_ids] =...
            compute_spatial_sensitivity(datarun, cone_finding_cell_spec, spat_sens_params);
        
        
        
        
        

        % IDENTIFY CONE LOCATIONS AND COLORS

        if 1

            switch combine_rfs

                case 'cluster'
                    
                    % FIND CONES LOCALLY WITHIN EACH RF
                    
                    % get local maxima from each RF
                    [initial_cone_centers_all,cone_sources] = find_cones_in_each_rf(...
                        datarun,spatial_cell_ids,all_sig_stixels,...
                        'strength',spat_sens_params.strength,'selection_params',spat_sens_params.selection_params,...
                        'filter',spat_sens_params.filter,'verbose',1);

                    % cluster local maxima
                    [clusters,cluster_centers,cluster_sources,max_dists] = ...
                        cluster_cone_centers(datarun,initial_cone_centers_all,cone_sources,cluster_rf_params);
                    
                    % fit locations, colors
                    cone_rf_params.cone_remap = remap;
                    cone_rf_params.cone_weights_figure = 20;
                    tic;
                    [cone_rgb, cone_spatial_profiles, cone_ideal_shapes, cone_rfs, cone_ids] = ...
                        summarize_cone_rfs2(datarun, cluster_centers, ...
                        cluster_sources,cone_kernel_params, cone_rf_params);
                    toc;

                    % get list of fitted cone centers
                    cone_centers = zeros(length(cone_ideal_shapes),2);
                    for cc = 1:length(cone_ideal_shapes)
                        cone_centers(cc,1:2) = cone_ideal_shapes{cc}.center;
                    end
                    
                    
                    
                    
                case 'spat sens'
                    % find local maxima in the spat. sens., then fit to each RF
                    
                    
                    % FIND LOCAL MAXIMA IN THE SPATIAL SENSITIVTY

                    % use the accumulated spatial sensitivity to find the number of cones and their approximate locations

                    switch cone_finding

                        case 'segment'  % identify pixels composing each cone using contiguity of significant pixels
                            cones_labeled = bwlabel(spatial_sensitivity,8);

                        case 'maxima'  % find local maxima, then segment

                            % find local maxima
                            local_maxima = find_cones_in_rf(spatial_sensitivity,'filter',[],...
                                struct('selection',struct('type','max','thresh',0)));

                            % segment
                            cones_labeled = bwlabel(local_maxima,8);

                        case 'manual'

                            % load user center points, identify which are sampled by which cell
                            [all_sig_stixels,cones_labeled,initial_cone_centers] = ...
                                use_cone_locations(datarun,cone_finding_cell_spec,cone_marks,'verbose',1);
                    end


                    % get number of cones
                    num_all_cones = max(max(cones_labeled));


                    % get cone center points
                    switch cone_finding

                        case {'segment','maxima'}
                            initial_cone_centers = zeros(num_all_cones,2);
                            for nn = 1:num_all_cones
                                [initial_cone_centers(nn,1),initial_cone_centers(nn,2)] = ait_centroid(cones_labeled==nn);
                            end
                    end



                    % GET PRECISE LOCATION AND COLOR PROFILE OF EACH CONE
                    % use individual STAs to refine the location of cones and extract their color profiles
                    
                    % fit BW version of initial cone locations, no penalty

                    switch which_cones
                        case 'all'  % use all cones
                            cone_analysis_roi = [];

                        case 'roi'  % create a ROI if it doesn't already exist
                            if ~exist('cone_analysis_roi','var')
                                figure;clf;cone_analysis_roi = roipoly(spatial_sensitivity);
                            end
                    end

                    % use roi
                    cone_rf_params.roi = cone_analysis_roi;

                    % pass in cone remapping function (used only if points are plotted)
                    cone_rf_params.cone_remap = remap;

                    % compute cone locations and colors
                    [cone_rgb, cone_spatial_profiles, cone_ideal_shapes, cone_rfs, cone_ids] =...
                        summarize_cone_rfs(datarun, spatial_cell_ids, initial_cone_centers,...
                        cones_labeled, all_sig_stixels, cone_kernel_params, cone_rf_params);

                    % get list of fitted cone centers
                    cone_centers = zeros(length(cone_ideal_shapes),2);
                    for cc = 1:length(cone_ideal_shapes)
                        cone_centers(cc,1:2) = cone_ideal_shapes{cc}.center;
                    end



            end
        end % If 1 (IDENTIFY CONE LOCATIONS AND COLORS)


        
        
        % GET ROI FOR CONE MOSAIC
        % estimate density, and only include cones in which local density is close to average density
        
        if 1
            
            % parameters
            num_bins = 20;
            bin_width = 1;

            % get local density
            [drp,bin_centers,drp_extras] = density_recovery_profile(cone_centers,num_bins,bin_width);

            % identify ROI
            cone_roi = identify_cone_mosaic_roi(cone_centers,drp_extras.density, drp_extras.eff_rad);
            
        end
            


        % CLASSIFY CONES
        % use cone color profiles to classify each as L, M, S, or unknown

        if 1
            % classify cones
            [cone_types,likelihood_ratios,classify_extras] = ...
                classify_cones(cone_rgb, cone_rgb_expected(datarun),classify_params);
        end



        % IDENTIFY CONE SAMPLING BY EACH RGC
        % explain each RF as a weighted sum of identified cones

        if 1
            % free up memory prior to this computation
            clear all_sig_stixels cone_ideal_shapes

            [cone_weights, regr_cell_ids,Wc]= extract_cone_weights(datarun, regression_cell_spec,...
                cone_spatial_profiles, cone_types, classify_extras.cone_rgb_observed);
            
            % clear RFs afterwards
            datarun.stas.rfs = cell(length(datarun.cell_ids),1);
        end



        % FIT CENTER AND SURROUND GAUSSIANS TO EACH RGC

        if 1
            clear cone_info
            cone_info.cone_weights = cone_weights;
            cone_info.cone_centers = cone_centers;
            cone_info.cell_ids = regr_cell_ids;

            rf_cone_fits = fit_cone_rfs(datarun,regression_cell_spec,...
                'cone_info',cone_info,'foa_profile',[],'fit_radius',150,'verbose',1,'foa_2d',4);

        end







%%      SAVE RESULTS


        % SAVE STUFF TO DISK
        %
        % in mat files and text files
        
        if save_results

            
            % EXPORT CONE SAMPLING DATA
            % write text files to disk
            % after writing to disk, load the data into datarun!
        
            if export_rgc_roi
                % create a ROI if it doesn't already exist
                if ~exist('rgc_roi','var') || isempty(rgc_roi)
                    % plot cone mosaic in a new figure
                    plot_cone_mosaic(datarun, 'fig_or_axes',0);
                    rgc_roi = roipoly;
                end
            else
                % export all cells
                rgc_roi = [];
            end

            % load info into a struct
            all_data.cone_weights = cone_weights;
            all_data.cone_types = cone_types;
            all_data.cone_centers = cone_centers;
            all_data.cone_rgb = cone_rgb;
            all_data.rf_cone_fits = rf_cone_fits;
            
            if isfield(classify_extras,'types_kmeans')
                all_data.cone_types_kmeans = classify_extras.types_kmeans;
                all_data.cone_types_em = classify_extras.types_em;
                all_data.cone_likelihoods = likelihood_ratios;
            else
                all_data.cone_types_kmeans = cone_types;
                all_data.cone_types_em = cone_types;
                all_data.cone_likelihoods = ones(length(cone_types),1);
            end
                

            % save most info to text files
            export_single_cone_data(datarun,regression_cell_spec,all_data,save_file_path,'rgc_roi',rgc_roi,'cone_roi',cone_roi)
            clear all_data

            % load data into datarun (needed for plotting stuff below)
            datarun = import_single_cone_data(datarun,save_file_path);
            


            % MAKE REPORT FIGURES
            % for quick visualization of the dataset

            % choose figures
            fig_mosaics = 110;
            fig_profiles = 111;
            fig_cone_mosaic = 112;
            fig_cone_classification = 113;
            fig_cone_spacing = 114;

            
            % spatial sensitivity, cell sampling, classification
            
            % set up plot axes
            figure(fig_mosaics);clf;
            pa.cone_field = subplot('Position',[.1 .7 .8 .25]);
            pa.grid = subplot_axes(fig_mosaics,[0 0 1 0.65],.05,.05,3,2);

            % plot things there

            % image of spatial sensitivity
            axes(pa.cone_field);imagesc(matrix_scaled_up(spatial_sensitivity,4));axis image;
            title(gca,printable_name(datarun))
            set(gca,'XTick',[],'YTick',[])

            % compute mosaic properties
            datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
            
            % mosaic sampling of each cell type
            for cc=1:5
                plot_cell_sampling(datarun, {cc},'fig_or_axes',pa.grid{cc},...
                    'cone_size',1.5,'line_width',[.05 .3]);
                set(gca,'XTick',[],'YTick',[])
                title(gca,datarun.cell_types{cc}.name)
            end

            % cone classification
            plot_cone_classification(datarun,'foa_hist',pa.grid{6},'dot_size',3);
            xlabel('');ylabel('')

            % print
            print(fig_mosaics,[save_file_path 'report'],'-dpdf')
            
            
            
            % average profiles
            
            % plot
            plot_average_rf_cone_profiles(datarun,1:5,'by_cone_type',1,...
                'figure',fig_profiles,'bin_edges', 0:2.5:100,'average_y','median','normalize','fit')
            % print
            print(fig_profiles,[save_file_path 'profiles'],'-dpdf')


            % cone mosaics and classification figures
            
            % plot 
            plot_cone_summary(datarun,'fig_drp',fig_cone_mosaic,'fig_class',fig_cone_classification,'cone_roi',cone_roi);
            % print
            print(fig_cone_mosaic,[save_file_path 'mosaic'],'-dpdf')
            print(fig_cone_classification,[save_file_path 'classification'],'-dpdf')
            

            
            % SAVE A FEW VARIABLES
            
            % spatial sensitivity spatial, if applicable
            save([save_file_path 'spatial_sensitivity'],'spatial_sensitivity')

            % save Wc
            save([save_file_path 'Wc'],'Wc')
            clear Wc

            % save cone RFs to matlab file
            avg_cone = compute_average_cone(cone_rfs);
            light_path_info = [datarun.piece.rig '-' datarun.piece.optical_path_direction];
            save([save_file_path 'cone_rfs'],'cone_rfs','avg_cone','light_path_info')
            clear cone_rfs
            
            % cluster stuff
            save([save_file_path 'cluster_stuff'],'initial_cone_centers_all','cone_sources',...
                'clusters','cluster_centers','cluster_sources','max_dists')
            
            
            
            
            % save analysis parameters
            all_params = struct;
            all_params.marks_params = marks_params;
            all_params.cone_finding_cell_spec = cone_finding_cell_spec;
            all_params.regression_cell_spec = regression_cell_spec;
            all_params.which_cells = which_cells;
            all_params.which_cones = which_cones;
            all_params.blur_rfs = blur_rfs;
            all_params.blur_radius = blur_radius;
            all_params.spat_sens_params = spat_sens_params;
            all_params.combine_rfs = combine_rfs;
            all_params.cluster_rf_params = cluster_rf_params;
            all_params.cone_finding = cone_finding;
            all_params.cone_rf_params = cone_rf_params;
            all_params.cone_kernel_params = cone_kernel_params;
            all_params.classify_params = classify_params;
            save([save_file_path 'parameters'],'all_params')

            
            % plot NND and latticing diagram
            figure(fig_cone_spacing);clf
            
            subplot(221);plot(mod(datarun.cones.centers(:,1),1),mod(datarun.cones.centers(:,2),1),'.');
            axis equal
            subplot(222);hist([mod(datarun.cones.centers(:,1),1) mod(datarun.cones.centers(:,2),1)]); 
            nnd = ipdm(datarun.cones.centers,'subset','nearest');
            nd = reshape(nnd,[],1);
            subplot(223);hist(nd(nd>0),0:.1:10);xlim([0 10])
            print(fig_cone_spacing,[save_file_path 'nnd_lattice'],'-dpdf')
            
            

        end
        
        



        % PLOT CONE SAMPLING BY INDIVIDUAL CELLS
        % spider plot

        if 0
            plot_cell_sampling(datarun, {4},'fig_or_axes',110)

            %plot_cell_sampling(datarun, datarun.cell_ids(get_cell_indices_roi(datarun,{2},output_roi)),...
            % cone_weights, regr_cell_ids, cone_centers, cone_types,struct('fig_or_axes',13))
        end



        % COMPARE MEASURED RFs TO RECONSTRUCTED CONE RFs
        % generate serial plots

        if 0

        end



%%



    catch
        
        % get the error
        temp=lasterror;
        
        % display the error and some eye-catching text
        fprintf('\n\n\n******* ERROR in %s *******\n\n',short_name)
        disp(temp.message)
        disp(temp.identifier)
        for ss =1:length(temp.stack)
            disp(temp.stack(ss))
        end
        fprintf('\n\n\n')

        % give an error so the script stops
        if ~catch_errors;
            rethrow(lasterror)
        end
        
    end
end








% PLOT CELLS WITH COLORED CONES

% if 0
%
%     % go through list of cells
%     for cc = 1:length(plot_cell_nums)
%
%         clear rf
%
%         % construct picture
%         rf(:,:,1) = reshape(cone_profile(:,L_cones)*cone_weights(L_cones,cc),...
%               datarun.stimulus.field_height,datarun.stimulus.field_width);
%         rf(:,:,2) = reshape(cone_profile(:,M_cones)*cone_weights(M_cones,cc),...
%               datarun.stimulus.field_height,datarun.stimulus.field_width);
%         rf(:,:,3) = reshape(cone_profile(:,S_cones)*cone_weights(S_cones,cc),...
%               datarun.stimulus.field_height,datarun.stimulus.field_width);
%
%         % plot original cell in first panel, reconstructed in second
%
%         % set up figure
%         figure(111);clf
%
%         % plot original
%         subplot(121);
%         imagesc(norm_image(datarun.stas.rfs{plot_cell_nums(cc)}))
%         axis image
%         title(num2str(datarun.cell_ids(plot_cell_nums(cc))))
%
%         % plot reconstructed
%         subplot(122);
%         imagesc(norm_image(rf));axis image;hold on
%
%
%         % add lines showing which cones are sampled
%
%         % get cell center point
%         rf_ctr = datarun.stas.rf_coms{plot_cell_nums(cc)};
%
%         % get list of strong cones
%         strong_cones = find(abs(cone_weights(:,cc)) > 6*robust_std(cone_weights(:,cc)));
%
%         % draw line from center point to each strong cone
%         for ss = 1:length(strong_cones)
%
%             % get cone center point
%             cone_ctr = cone_centers(strong_cones(ss),:);
%
%             % plot line
%             plot([rf_ctr(1) cone_ctr(1)],[rf_ctr(2) cone_ctr(2)],'k')
%         end
%
%
%         % link axes
%         linkaxes([subplot(121) subplot(122)])
%
%         % zoom in on the cell
%         zoom_rad = 35;
%         set(gca,'XLim',[rf_ctr(1)-zoom_rad rf_ctr(1)+zoom_rad],'YLim',[rf_ctr(2)-zoom_rad rf_ctr(2)+zoom_rad ])
%
%         pause
%
%     end
%
% end


% remove java path
%javarmpath /snle/lab/Applications/Vision.app/Contents/Resources/Java/Vision.jar
%clear java

