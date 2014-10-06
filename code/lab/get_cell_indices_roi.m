function [cell_indices,cells_kept] = get_cell_indices_roi(datarun, cell_spec, roi, varargin)
% get_cell_indices_roi     Return only cells within a particular ROI
%
% usage:  cell_indices = get_cell_indices_roi(datarun, cell_spec, roi, <params>)
%
% arguments:  datarun - datarun struct
%           cell_spec - cell specifiation (see get_cell_indices for options)
%                 roi - region of interest.  only return cells with a center point inside this region.
%                           if logical, must be same size as stas, and is interpretted as logical map
%                           if struct, is interpretted as polygon in standard format (see PolygonClipper documentation)
%                           if empty, all cells are returned
%              params - struct or list of optional parameters (see below)
%
% outputs:
%        cell_indices - numbers of cells with center points inside the ROI
%          cells_kept - logical vector indicating which of the possible cells were kept
%
%
% optional fields in params, their default values, and what they specify:
%
% cell_locations     	'center'      	how to locate each cell
%                                       'center' - use the center point, as returned by function 'rf_centers'
%                                       'marks' - use the marks.  if any mark falls within the ROI, keep the cell
%
%
% 2009-08  gauthier
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('cell_locations', 'center', @(x)any(strcmpi(x,{'center','marks'})));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;






% CHECK FOR ERRORS IN INPUT
switch class(roi)
    case 'logical'
        % logical matrix must match the size of the STAs
        % first ensure it is a matrix
        if length(size(roi)) ~= 2 || ...
                ... then ensure it is the same size
                ~all(size(roi) == [datarun.stimulus.field_height datarun.stimulus.field_width])
            error('Size of ROI does not match size of STAs')
        end

    case 'struct'
        % ensure it is a polygon

        % THIS FUNCTIONALITY SHOULD BE CREATED IN A SEPARATE FUNCTION, SOMETHING LIKE 'assert_polygon'
end





% LOAD LOCATIONS OF ALL CELLS IN THE STYLE OF MARKS (I.E. STRUCT OF LOGICAL MATRICES)

% get cell indices for all cells
all_cell_indices = get_cell_indices(datarun,cell_spec);

% if empty, return all cells
if isempty(roi)
    cell_indices = all_cell_indices;
    cells_kept = true(length(cell_indices),1);
    return
end

% initialize cell location storage
cell_locations = cell(length(all_cell_indices),1);

% fill in cell_locations
for cc = 1:length(all_cell_indices)

    % initialize as empty
    cell_locations{cc} = false(datarun.stimulus.field_height,datarun.stimulus.field_width);

    % load up location(s) for this cell
    switch params.cell_locations
        case 'center'
            
            % get cell center (will return [-1 -1] if center doesn't exist)
            center = rf_center(datarun,datarun.cell_ids(all_cell_indices(cc)));
            % round to be integers
            center = round(center);
            % if center exists...
            if ~isempty(center) && all(center>0)
                % enter it
                cell_locations{cc}(center(1,2),center(1,1)) = true;
            end

            
        case 'marks'

            % get marks, if they exist
            if ~isempty(datarun.stas.marks{all_cell_indices(cc)})
                cell_locations{cc} = datarun.stas.marks{all_cell_indices(cc)};
            end

    end
end





% FIND WHICH CELLS HAVE LOCATIONS IN THE ROI


% start list of cells to be kept
cells_kept = false(length(all_cell_indices),1);

% check each cell
for cc = 1:length(all_cell_indices)
    switch class(roi)
        case 'logical'

            % keep if any intersection
            if any(any(cell_locations{cc} & roi))
                cells_kept(cc) = true;
            end

        case 'struct'   % interpret as polygon

            %  NOT SUPPORTED YET
            error('ROI from polygon not supported yet.')

    end
end



% note cells to keep
cell_indices = all_cell_indices(cells_kept);







% % should be replaced with "center = rf_center(datarun,cell_id,center_type)"
% for cc = 1:length(all_cell_nums)
%     switch params.center_type
%         case 'rf_com'
%             % get center point (might be empty)
%             cp = datarun.stas.rf_coms{all_cell_nums(cc)};
%
%             % if non empty, keep it
%             if ~isempty(cp)
%                 center_points(cc,1:2) = cp;
%             else
%                 % if empty, set this value to NaN NaN
%                 center_points(cc,1:2) = [NaN NaN];
%             end
%     end
% end
%
%
% % remove cells from list which had no location
%
% % find those cells
% cells_to_check = find(~isnan(center_points(:,1)));
%
% % exclude them
% cell_nums = all_cell_nums(cells_to_check);
% center_points = center_points(cells_to_check,1:2);


%
% % alternative one-liner
% cells_to_keep = find(diag(roi(round(center_points(:,2)),round(center_points(:,1)))));
%
%
% if 0
%     % get value of ROI at cell center point (rounded)
%     roi_val = roi(round(center_points(cc,2)),round(center_points(cc,1)));
%
%     % if value is 1
%     if roi_val
%         cells_to_keep = [cells_to_keep cc];
%     end
%
% end


%
% % return cells to keep
% cell_indices = cell_nums(cells_to_keep);
%
%
% % return logical vector of which cells were kept
% if nargout > 1
%     % initialize vector
%     cells_kept = false(length(all_cell_nums),1);
%     % go through list of all cell nums
%     for cc = 1:length(all_cell_nums)
%         % if it was in the list of returned cell numbers...
%         if ~isempty(intersect(cell_indices,all_cell_nums(cc)))
%             % set the value to true
%             cells_kept(cc) = true;
%         end
%     end
% end



