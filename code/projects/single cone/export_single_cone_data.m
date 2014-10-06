function export_single_cone_data(datarun,cell_spec,all_data,file_path,varargin)
% EXPORT_SINGLE_CONE_DATA     Create text files containing single cone sampling data
%
% usage:  export_single_cone_data(datarun,cell_spec,cone_weights,...
%             cone_types,cone_centers,cone_rgb,rf_cone_fits,file_path,varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells to write to disk
%            all_data - struct with fields 'cone_weights','cone_types','cone_centers','cone_rgb','rf_cone_fits'
%                           if empty, this info is taken from the appropriate fields in datarun
%           file_path - prefix of where to save
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional fields in params, their default values, and what they specify:
%
% verbose           true        show what information is being saved
% rgc_roi           []          ROI for RGCs.  format like an argument to get_cell_indices_roi.
% cone_roi          []          Nx1 boolean vector.  empty is equivalent to all true
%
%
%
%
% these are the columns present in each file
% (items with *** are not implemented yet)
% 
% 
% rgcs.txt
% 
% cell id
% cell type (0 = none, 1 = on-parasol, 2 = off-parasol, 3 = on-midget, 4 = off-midget, 5 = sbc, 6+ = other)
% ROI (1 = in ROI, 0 = not)
% cone COM x  (if COM doesn't exist, [-1 -1] is used) 
% cone COM y
% RGC RF gaussian fit (based on single cone RF) mean, x
% RGC RF gaussian fit (based on single cone RF) mean, y
% RGC RF gaussian fit (based on single cone RF) center radius
% RGC RF gaussian fit (based on single cone RF) center scale
% RGC RF gaussian fit (based on single cone RF) surround radius
% RGC RF gaussian fit (based on single cone RF) surround scale
% radius in which the fit was computed (units = stixels)
% error of the fit
% ***contamination
% ***num-spikes
% 
% 
% cones.txt
% 
% cone ID
% x location
% y location
% identity (0 = unsure, 1 = L, 2 = M, 3 = S)
% cone R sensitivity (units proportional to strength)
% cone G sensitivity
% cone B sensitivity
% ROI (1 = in ROI, 0 = not)
% identity derived from EM cone classification
% likelihood ratio from EM cone classification
% identify derived from kmeans cone classification
% ***sd (sd of Gaussian fit probably constant but just in case)
% ***identity-confidence (for starters, let's just put in the L/M likelihood ratio here, and 0.0 for the S cones)
% 
% 
% connections.txt
% 
% weight matrix of connection strengths
%    rows = cones, columns = cells
%    
%
%
%
% examples:
%
%
%   read from datarun:
%
%     export_single_cone_data(datarun,{1,2,3,4,5},[],[single_cone_path datarun.names.short_name '/']);
%
%     
%   read from all_data:
%
%     % load info into a struct
%     all_data.cone_weights = cone_weights;
%     all_data.cone_types = cone_types;
%     all_data.cone_centers = cone_centers;
%     all_data.cone_rgb = cone_rgb;
%     all_data.rf_cone_fits = rf_cone_fits;
% 
%     % save most info to text files
%     export_single_cone_data(datarun,regression_cell_spec,all_data,save_file_path)
%    
% gauthier, 2008-09
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', true);
p.addParamValue('rgc_roi', []);
p.addParamValue('cone_roi', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;






% SPECIFY WHAT DATA TO WRITE OUT

% make list of which pieces of information to include in each file
% first item is the file name, subsequent items are the columns to include
rgcs_file = {[file_path 'rgcs.txt'],'cell_ids','cell_type','rgc_roi','rgc_com','rf_cone_fits'};
cones_file = {[file_path 'cones.txt'],'cone_ids','cone_centers','cone_types','cone_rgb','cone_roi','cone_types_em','cone_likelihood','cone_types_kmeans'};
connections_file = {[file_path 'connections.txt'],'cone_weights'};

% put all files into a list of files
file_list = {cones_file,rgcs_file,connections_file};



% PUT CONE INFORMATION INTO A STRUCT (all_data)

% if empty to begin with, fill in with single cone info from datarun
if isempty(all_data)
    
    % save information from all the cones...
    all_data.cone_types = datarun.cones.types;
    all_data.cone_centers = datarun.cones.centers;
    all_data.cone_rgb = datarun.cones.rgb;
    all_data.cone_types_kmeans = datarun.cones.types_kmeans;
    all_data.cone_types_em = datarun.cones.types_em;
    all_data.cone_likelihoods = datarun.cones.likelihoods;
    
    % ...but only the specified RGCs
    % get cell_indices for cells which are to be returned
    cell_indices = get_cell_indices(datarun,cell_spec);
    % only save fits and weights from those cells
    all_data.rf_cone_fits = {datarun.cones.rf_fits{cell_indices}};
    all_data.cone_weights = datarun.cones.weights(:,cell_indices);
end

% store other input information in the same struct
all_data.datarun = datarun;
all_data.cell_spec = cell_spec;
all_data.params = params;



% COLLECT INFORMATION FOR EACH FILE AND WRITE IT OUT


% initialize struct of file names that were written out
filenames = cell(0);

% get numbers and write out each file
for ff = 1:length(file_list)
    
    % get list for this file
    this_file = file_list{ff};
    
    % initialize matrix to write out
    out_matrix = [];
    
    % go through each piece of information and include it in the matrix
    for gg = 2:length(this_file)
        if params.verbose
            fprintf('getting %s...\n',this_file{gg})
        end
        temp_content = get_content(all_data,this_file{gg});
        % if the field is empty fill it with zeros -- this was introduced by gdf to make compatable with new cone finding - 
        % i.e. no EM based cone classification
        if isempty(temp_content)
            element_num = size(out_matrix,1);
            temp_content(1:element_num,1) = 0;
        end
        out_matrix = [out_matrix temp_content];
    end
    
    % write out file, tab delimited
    dlmwrite(this_file{1}, out_matrix, '\t')

    % note the file name
    filenames{ff} = this_file{1};
    
    % display the file name
    if params.verbose
        fprintf('     [%d x %d] matrix written to file ''%s''...\n\n',...
            size(out_matrix,1),size(out_matrix,2),this_file{1})
    end
    
end




% FUNCTION SPECIFYING WHAT CONTENT TO GET

function content = get_content(ad,request)
% ad - struct containing all data
% request - string specifying what information is desired

switch request
    case 'cell_ids'
        % return cell_ids in a column vector
        content = reshape(ad.datarun.cell_ids(get_cell_indices(ad.datarun,ad.cell_spec)),[],1);
        
    case 'rgc_roi'
        % return 1 for RGCs in the ROI, 0 for not
        
        % get_cell_indices_roi does all the work
        [cell_numbers,cells_kept] = get_cell_indices_roi(ad.datarun,ad.cell_spec,ad.params.rgc_roi);
        
        % return just the logical vector of which cells were in the ROI
        content = cells_kept;

    case 'cell_type'
        % return list of cell types for each cell in a column vector
        cell_ids = ad.datarun.cell_ids(get_cell_indices(ad.datarun,ad.cell_spec));
        content = reshape(find_cell_types(ad.datarun,cell_ids),[],1);

    case 'rgc_com'
        % return list of COM x and y coord for each RGC
        cell_indices = get_cell_indices(ad.datarun,ad.cell_spec);
        
        % initialize variable
        coms = zeros(length(cell_indices),2);
        
        % get COM of each cell
        for cc = 1:length(cell_indices)
            
            % get COM
            com_temp = ad.datarun.stas.rf_coms{cell_indices(cc)};
            
            % if COM is empty, return [-1 -1]
            if isempty(com_temp)
                coms(cc,1:2) = [-1 -1];
            else
                
                % if it does exist, ensure it is positive, or equal to [-1 -1]
                if (com_temp(1) > 0 && com_temp(2) > 0) || all(com_temp == -1)
                    coms(cc,1:2) = com_temp;
                else
                    
                    % otherwise give an error
                    error('COM for cell id %d is non-positive (COM is [%0.2f %0.2f]).',...
                        ad.datarun.cell_ids(cell_indices(cc)),com_temp(1),com_temp(2))
                end
            end
        end
        
        % return the COMs
        content = coms;

    case 'cone_ids'
        % return list of cone ids
        content = [1:length(ad.cone_types)]';

    case 'cone_centers'
        % return cone center locations
        content = ad.cone_centers;

    case 'cone_types'
        % return cone types
        content = str2num(strrep(strrep(strrep(strrep(ad.cone_types','L','1'),'M','2'),'S','3'),'U','0')');

    case 'cone_types_em'
        % return cone types
        content = str2num(strrep(strrep(strrep(strrep(ad.cone_types_em','L','1'),'M','2'),'S','3'),'U','0')');
 
    case 'cone_types_kmeans'
        % return cone types
        content = str2num(strrep(strrep(strrep(strrep(ad.cone_types_kmeans','L','1'),'M','2'),'S','3'),'U','0')');
 
    case 'cone_rgb'
        % return rgb triplet of each cone
        content = ad.cone_rgb;
        
    case 'cone_roi'
        % return 1 for cones in the ROI, 0 for not
        
        % if empty, all cones are included
        if isempty(ad.params.cone_roi)
            cells_kept = true(size(ad.cone_centers,1),1);
        else
            % if non empty, just the roi
            cells_kept = ad.params.cone_roi;
        end
        
        % return the logical vector of which cones were in the ROI
        content = cells_kept;
        
    case 'cone_likelihood'
        % likelihood ratio for each cone
        
        content = ad.cone_likelihoods;

    case 'cone_weights'
        % return weight of each cone
        content = ad.cone_weights;

    case 'rf_cone_fits'
        % return all parameters of a DOG fit to each single cone RF
        
        % get number of fits
        num_fits = length(ad.rf_cone_fits);
        
        % initialize output matrix
        content = zeros(num_fits,8);
        
        for ff = 1:num_fits
            
            % get the fit
            fit = ad.rf_cone_fits{ff};
            
            % if it's not empty...
            if ~isempty(fit)
                % return the following parameters
                content(ff,:) = [...
                    fit.center...
                    fit.center_radius...
                    fit.center_scale...
                    fit.surround_radius...
                    fit.surround_scale...
                    fit.fit_radius...
                    fit.error];
            end
        end
        
        
        
    case 'fit_center_temp'
        % return center of the vision fit, unfiltered
        % THIS IS A POOR CHOICE OF CENTER POINT, ONLY TO BE USED TEMPORARILY

        % get list of cell numbers
        cell_indices = get_cell_indices(ad.datarun,ad.cell_spec);

        % initialize variable
        fit_centers = zeros(length(cell_indices),2);

        % get fit center of each cell
        for cc = 1:length(cell_indices)
            fit_centers(cc,1:2) = rf_center(ad.datarun,ad.datarun.cell_ids(cell_indices(cc)),'vision');
        end
        
        content = fit_centers;

    case 'fit_radius_temp'
        % return radius of the vision fit, unfiltered
        % THIS IS A POOR CHOICE OF RADIUS, ONLY TO BE USED TEMPORARILY

        % get list of cell numbers
        cell_indices = get_cell_indices(ad.datarun,ad.cell_spec);

        % initialize variable
        fit_radii = zeros(length(cell_indices),1);

        % get fit center of each cell
        for cc = 1:length(cell_indices)
            fit_radii(cc) = sqrt(prod(ad.datarun.vision.sta_fits{cell_indices(cc)}.sd));
        end
        
        content = fit_radii;

    case 'name'
        % return 
        content = [];

    otherwise
        error('don''t know how to get information type ''%s''.',request)
end


% ensure it's double
content = double(content);


% contamination
% num-spikes
%
% File: data000-cones
% cone-ID
% x
% y
% sd (sd of Gaussian fit probably constant but just in case)
% identity (L,M,S)
% identity-confidence (for starters, let's just put in the L/M likelihood ratio here, and 0.0 for the S cones)
% SNR (of the signal from the cell it was identified)
%
% File: data000-connections
% rgc-ID
% cone-ID
% weight (expressed in units relative to other cones providing input to the same rgc, i.e. like the STA )
