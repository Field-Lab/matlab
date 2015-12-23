function profiles = get_nearest_neighbor_profile(datarun, cell_spec, varargin)
% get_nearest_neighbor_profile     identify average RF profiles along the line
%                                  connecting nearest neighboring RF centers
%
% usage:  profiles = get_nearest_neighbor_profile(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:   profiles - 2xNxC matrix, N = # samples, C = # colors in the RF
%                           row 1 is the reference cell, row 2 is the nearest neighbor
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           true    display status
% extension         3       how much to extend the segment beyond cell centers 
% samples           100     how many points to sample along the line connecting neighboring cell centers
%
%
%
% 2009-04 gauthier
%


% TO ADD
%
% * which center point to use
% * which RF to use
%




% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('extension', 3);
p.addParamValue('samples', 100);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% LOAD CENTER POINTS

% get center points of cells
[centers,cell_ids] = rf_centers(datarun,cell_spec);

% get cell indices of the cells with center points
cell_indices = get_cell_indices(datarun,cell_ids);



% IDENTIFY NEAREST NEIGHBORS

nbrmat = ipdm(centers,'result','structure','subset','nearestNeighbor');



% COLLECT RF SLICES ALONG NEAREST NEIGHBOR LINES


% show output
if params.verbose
    fprintf('\nExtracting slices for %d cells...',length(cell_indices));
    start_time = clock; % note when it started
end

% initialize variables storing slices
ref_slices = [];
nbr_slices = [];

% loop through cells
for cc = 1:length(cell_indices)
    
    % get cell index and id
    cell_index = cell_indices(cc);
    cell_id = datarun.cell_ids(cell_index);
    
    % identify nearest neighbor
    cell_index_nbr = cell_indices(nbrmat.columnindex(cc));
    cell_id_nbr = datarun.cell_ids(cell_index_nbr);
    
    % get RFs
    rf = get_rf(datarun,cell_id);
    rf_nbr = get_rf(datarun,cell_id_nbr);
    
    % verify RFs exist
    if isempty(rf) || isempty(rf_nbr)
        continue
    end
    
    % define segment connecting centers
    segment = centers([cc nbrmat.columnindex(cc)],:);
    
    % skip if duplicates -- THIS SHOULD BE DONE AHEAD OF TIME, NOT HERE!!!
    %if norm(segment(1,:)-segment(2,:))<1; continue; end
    
    % get slices
    ref_slice = matrix_slice(rf,segment,'extension',params.extension,'samples',params.samples);
    nbr_slice = matrix_slice(rf_nbr,segment,'extension',params.extension,'samples',params.samples);
    
    % zero out points beyond the edge so they don't interfere with computing the mean later
    ref_slice(isnan(ref_slice)) = 0;
    nbr_slice(isnan(nbr_slice)) = 0;
    
    % save slices
    ref_slices = [ref_slices; permute(ref_slice,[3 1 2])];
    nbr_slices = [nbr_slices; permute(nbr_slice,[3 1 2])];
    
    % debugging option: plot the slices
    %figure(6);clf;plot([ref_slice nbr_slice]);pause
end


% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end



 % AVERAGE THE SLICES TOGETHER
 
 % stupid matlab requires a loop for these assignments to work consistently!
 
%  profiles(1,:,:) = mean(ref_slices,1);
%  profiles(2,:,:) = mean(nbr_slices,1);
 
for cc=1:size(ref_slices,3)
    profiles(1,:,cc) = mean(ref_slices(:,:,cc),1);
    profiles(2,:,cc) = mean(nbr_slices(:,:,cc),1);
end
