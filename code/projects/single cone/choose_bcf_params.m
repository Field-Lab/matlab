
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   BAYESIAN CONE FINDING STEP 2 OF 5    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose parameters for the analysis




% GENERAL PARAMETERS

% initialize struct
bcf_params = struct;

% size of subsets
bcf_params.padding_x = 5;
bcf_params.padding_y = 5;
bcf_params.roi_x_size = 2*bcf_params.padding_x + 10;
bcf_params.roi_y_size = 2*bcf_params.padding_y + 10;

% don't recompute if not necessary
bcf_params.new_W = 0;
bcf_params.new_STAs = 0;

% cell spec
%cone_finding_cell_spec = 'all';
bcf_params.cone_finding_cell_spec = {1,2,3,4,5,6};

% which cell types to identify the sampling of
bcf_params.regression_cell_spec = {1,2,3,4,5,6};
%bcf_params.regression_cell_spec = 'all';
%bcf_params.regression_cell_spec = datarun.cell_types{4}.cell_ids(1:5);

% figure
bcf_params.plot_fig = 20;
bcf_params.cones_fig = 21;
bcf_params.dll_fig = 22;

% iterations per patch
bcf_params.num_iter = 100 ;

% radius of relevance
% sets which stixels (around the marks) are considered relevant in each STA
bcf_params.rel_radius = 4;

% kernel RGB
bcf_params.kernel_colors = cone_rgb_expected(datarun);

% kernel plot colors
bcf_params.kernel_plot_colors = ('rgb')';
% TO CHANGE THIS: add a line to the dataset-specific parameters below

% start time of the datarun
bcf_params.start_time = start_time;


% DATASET-SPECIFIC PARAMETERS

switch datarun.names.nickname

    case 'blueberry'

        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % prior cone density
        bcf_params.q = 0.9; % blueberry

        % kernel spacing (in pixels) and radius
        %bcf_params.kernel_spacing = 1/2; bcf_params.kernel_radii = 0.75*[1 1 1];
        %bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.70*[1 1 1];
        bcf_params.kernel_spacing = 1/4;  bcf_params.kernel_radii = [.63 .63 .70];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [2.3 2.5]; bcf_params.S_S = [3.8 4];

        % central square
        %bcf_params.relevant_region_x = [82 102]; bcf_params.relevant_region_y = [103 123];
        % whole thing
        bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];

        

    case 'peach'

        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % prior cone density
        bcf_params.q = 0.05;

        % size, spacing of candidate cones
        bcf_params.kernel_spacing = 1/4;
        bcf_params.kernel_radii = 0.65*[1 1 1];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [1.8 2.0]; bcf_params.S_S = [3.8 4];

        % identify relevant rectangle
        switch 1
            case 1  % whole thing
                bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];
            case 2  % central square
                bcf_params.relevant_region_x = [120 193]; bcf_params.relevant_region_y = [127 201];
            case 3  % smaller central region
                bcf_params.relevant_region_x = [120 145]; bcf_params.relevant_region_y = [127 140];
            case 4  % ???
                bcf_params.relevant_region_x = [60 256]; bcf_params.relevant_region_y = [72 267];
        end

        
        
    case 'plantain'

        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % cone density prior
        bcf_params.q = 0.05;

        % kernel spacing (in pixels) and radius
        bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [1.8 2]; bcf_params.S_S = [3.8 4];

        % central square
        %bcf_params.relevant_region_x = [75 210]; bcf_params.relevant_region_y = [160 220];
        % whole thing
        bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];

    case 'plantain-1'
        
        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % cone density prior
        bcf_params.q = 0.05;

        % kernel spacing (in pixels) and radius
        bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [1.8 2]; bcf_params.S_S = [3.8 4];

        % central square
        %bcf_params.relevant_region_x = [75 210]; bcf_params.relevant_region_y = [160 220];
        % whole thing
        bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];
        
        
    case 'plantain_last_hour'

        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % cone density prior
        bcf_params.q = 0.05;

        % kernel spacing (in pixels) and radius
        bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [1.8 2]; bcf_params.S_S = [3.8 4];

        % central square
        %bcf_params.relevant_region_x = [75 210]; bcf_params.relevant_region_y = [160 220];
        % whole thing
        bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];

        
        
    case 'apricot'

        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % cone density prior
        bcf_params.q = 0.05;

        % kernel spacing (in pixels) and radius
        bcf_params.kernel_spacing = 1/5; bcf_params.kernel_radii = 0.65*[1 1 1];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [2.3 2.5]; bcf_params.S_S = [3.8 4];
        %bcf_params.LM_MM = [0.1 .5]; bcf_params.LM_S = LM_MM; bcf_params.S_S = [3.8 4];

        % central square
        %bcf_params.relevant_region_x = [160 480]; bcf_params.relevant_region_y = [50 270];
        % whole thing
        bcf_params.relevant_region_x = [1 640-10]; bcf_params.relevant_region_y = [1 320-10];

        
        
    case 'kiwi'

        % arbitrary scale factor
        bcf_params.magic_number = 1;
        
        % cone density prior
        bcf_params.q = 0.05;

        % kernel spacing (in pixels) and radius
        %bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.9*[1 1 1];
        bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = [0.9 0.9 0.7];

        % distance prior
        bcf_params.LM_MM = [2.8 3.1]; bcf_params.LM_S = [2.3 2.5]; bcf_params.S_S = [3.8 4];

        % whole thing
        bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];



    case 'grapes'

        % arbitrary scale factor
        bcf_params.magic_number = 1;
        
        % cone density prior
        bcf_params.q = 0.05;
        
        % kernel spacing (in pixels) and radius
        bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = [0.75 0.75 0.9];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [2.3 2.5]; bcf_params.S_S = [3.8 4];
        
        % whole thing
        bcf_params.relevant_region_x = [1 640-10]; bcf_params.relevant_region_y = [1 320-10];

        
    case 'butterfly'

        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % prior cone density
        bcf_params.q = 0.05; % blueberry

        % kernel spacing (in pixels) and radius
        %bcf_params.kernel_spacing = 1/2; bcf_params.kernel_radii = 0.75*[1 1 1];
        %bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.70*[1 1 1];
        bcf_params.kernel_spacing = 1/4;  bcf_params.kernel_radii = [.63 .63 .70];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [2.3 2.5]; bcf_params.S_S = [3.8 4];

        % whole thing
        bcf_params.relevant_region_x = [1 640-10]; bcf_params.relevant_region_y = [1 320-10];

        

    case 'apple-12'


        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % cone density prior
        bcf_params.q = 0.05;

        % kernel spacing (in pixels) and radius
        bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

        % distance prior
        bcf_params.C_C = [2.3 2.5]; 

        % central square
        % whole thing
        bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];

    

    case 'apple-13'

        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % cone density prior
        bcf_params.q = 0.05;

        % kernel spacing (in pixels) and radius
        bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [1.8 2]; bcf_params.S_S = [3.8 4];

        % central square
        %bcf_params.relevant_region_x = [75 210]; bcf_params.relevant_region_y = [160 220];
        % whole thing
        bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];


    case 'apple-15'


        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % cone density prior
        bcf_params.q = 0.05;

        % kernel spacing (in pixels) and radius
        bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

        % distance prior
        bcf_params.C_C = [2.3 2.5]; 

        % central square
        % whole thing
        bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];
        
        
    case 'plantain_15'

        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % cone density prior
        bcf_params.q = 0.05;

        % kernel spacing (in pixels) and radius
        bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [1.8 2]; bcf_params.S_S = [3.8 4];

        % central square
        %bcf_params.relevant_region_x = [75 210]; bcf_params.relevant_region_y = [160 220];
        % whole thing
        bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];

        
        
    case 'plantain_30'

        % arbitrary scale factor
        bcf_params.magic_number = 1;

        % cone density prior
        bcf_params.q = 0.05;

        % kernel spacing (in pixels) and radius
        bcf_params.kernel_spacing = 1/3; bcf_params.kernel_radii = 0.65*[1 1 1];

        % distance prior
        bcf_params.LM_MM = [2.3 2.5]; bcf_params.LM_S = [1.8 2]; bcf_params.S_S = [3.8 4];

        % central square
        %bcf_params.relevant_region_x = [75 210]; bcf_params.relevant_region_y = [160 220];
        % whole thing
        bcf_params.relevant_region_x = [1 320-10]; bcf_params.relevant_region_y = [1 320-10];


        
    otherwise
        error('datarun name %s not recognized!',datarun.names.nickname)
end




