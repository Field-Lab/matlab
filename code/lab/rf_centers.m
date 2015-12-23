function [centers,cell_ids] = rf_centers(datarun,cell_spec,center_type,coordinates)
% rf_center     return the center points of many RFs
%   
%   MODE 1 (nargout = 1): if no center point is found, [-1 -1] is returned as the center point
%
%   MODE 2 (nargout = 2): if no center point is found, the cell is excluded from the list
%                          	the second argument is the list of cell ids which are included
%
%
%    mode is set by the number of output arguments
%
%
% usage:  center = rf_centers(datarun,cell_spec,center_type)
%
%         [centers,cell_ids] = rf_centers(datarun,cell_spec,center_type)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells
%         center_type - specifies what kind of center point to return, see rf_center for options
%         coordinates - what coordinate system to use, see options in 'rf_center'
%
% outputs:    centers - Nx2 matrix of x,y coordinates
%            cell_ids - list of cell ids which are included (see above)
%
%
% gauthier  2009-02
%


% if no center type is specified, make empty
if ~exist('center_type','var')
    center_type = [];
end
% if no coordinates type is specified, make empty
if ~exist('coordinates','var')
    coordinates = [];
end

% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% initialize output matrix to be all -1
centers = -1*ones(length(cell_indices),2);

% loop through cells
for cc = 1:length(cell_indices)
    
    % get cell index and id
    cell_index = cell_indices(cc);
    cell_id = datarun.cell_ids(cell_index);
    
    % get the rf center
    temp = rf_center(datarun,cell_id,center_type,coordinates);
    
    % enter RF center, if is exists
    if ~isempty(temp) && ~all(temp==[1 1])
        centers(cc,:) = rf_center(datarun,cell_id,center_type,coordinates);
    end
        
end

% decide which mode
switch nargout
    case 1
        % nothing to be done
        
    case 2
        
        % only keep cells with a center point (i.e. their center is not [-1 -1])
        keep_indices = any(centers ~= -1,2);

        % get these cell ids
        cell_ids = datarun.cell_ids(get_cell_indices(datarun,cell_spec));
        cell_ids = cell_ids(keep_indices);

        % only keep these center points
        centers = centers(keep_indices,:);

    otherwise
        error('More than 3 argument not supported.  See help for info.')
end

