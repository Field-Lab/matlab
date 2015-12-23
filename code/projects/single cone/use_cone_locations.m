function [sig_stixels,cones_labeled,initial_cone_centers] = use_cone_locations(datarun,cell_spec,cone_marks,varargin)
% use_cone_locations     read in cone locations from an external source, such as user clicking
%
% usage:  [sig_stixels,cones_labeled,initial_cone_centers] = use_cone_locations(datarun,centers_file_path)
%
% arguments:          datarun - datarun struct
%                   cell_spec - which cells to use for cone finding
%                  cone_marks - cell struct, cone_marks{cell_index} is a Nx2 matrix of cone locations found in that cell
%
% outputs:    sig_stixels - 
%           cones_labeled - YxX matrix, equal to size of a single STA.  Mostly 0s, and at each cone location
%                               indicates the cone id
%    initial_cone_centers - 
%
%
% optional params, their default values, and what they specify:
%
% verbose           false           show output
% fig_or_axes       []              figure or axes to plot in. if 0, make new figure. if empty, don't plot
% thresh            4               threshold for considering a cone sampled by a cell
%                                       measured in units of SNR for that cell
%
% 2009-03 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('fig_or_axes', []);
p.addParamValue('thresh', 4);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;







% get cell indices for the given cell_spec
cell_indices = get_cell_indices(datarun,cell_spec);

% initialize variables
noise_sigmas = zeros(length(cell_indices),1);
initial_cone_centers = [];
sig_stixels = sparse(datarun.stimulus.field_height*datarun.stimulus.field_width,length(cell_indices));
cones_labeled = zeros(datarun.stimulus.field_height,datarun.stimulus.field_width);



% identify noise sigma of each RF

for cc=1:length(cell_indices)
    cell_index = cell_indices(cc);
    if isempty(datarun.stas.rfs{cell_index})
        error('RF of cell index %d not found',cell_index)
    end
    noise_sigmas(cell_index) = robust_std(reshape(datarun.stas.rfs{cell_index},[],1));
end


% build list of cone locations, and identify which cells sample from each cone


% show output
if params.verbose
    fprintf('\nIdentifying cell sampling for %d cells...',length(cone_marks));
    start_time = clock; % note when it started
end


% go through each cell
for mm=1:length(cone_marks);
    
    % go through each cone
    for cc=1:size(cone_marks{mm},1)

        % get the cell index and cell id of this cell
        cell_index = mm;
        cell_id = datarun.cell_ids(cell_index);

        % note cone location
        cone_loc = cone_marks{mm}(cc,:);

        % append to the list of cone locations
        initial_cone_centers=[initial_cone_centers;cone_loc];
        
        % store rounded version (just saves on typing below)
        round_loc = round(cone_loc);
        
        % check that no two cones occupy the same location
        if cones_labeled(round_loc(2),round_loc(1)) ~= 0
            % if they do, keep only the first
            fprintf(['\n\n\nWARNING!!! TWO CONE SEED LOCATIONS OCCUPY THE SAME PIXEL\n'...
                '  cell index %d, cone %d duplicates a previous cone location, and is being ignored\n\n\n'],...
                cell_index,cc)
            % delete this cone
            initial_cone_centers = initial_cone_centers(1:end-1,:);
            continue
        end
            
        % otherwise, enter the cone in cones_labeled
        cones_labeled(round_loc(2),round_loc(1)) = size(initial_cone_centers,1);

        
        % get value from each RF at this location
        
        % initialize list
        rf_vals = zeros(length(cell_indices),3);
        
        % fill in values
        for ci=1:length(cell_indices)
            % get RF
            rf = get_rf(datarun,datarun.cell_ids(cell_indices(ci)));
            % get normalized values
            rf_vals(ci,:) = rf(round_loc(2),round_loc(1),:)/noise_sigmas(cell_indices(ci));
        end
        
        
        % identify which cells sample signficantly from the cone
        
        % collapse across color by summing
        test_vals = sum(rf_vals,2);
        
        % find values that exceed threshold
        sig_cells = find(abs(test_vals) > params.thresh*robust_std(test_vals));
        
        % be sure to include the cell in which the user clicked on the cone
        %sig_cells = union(sig_cells,find(cell_indices==get_cell_indices(datarun,cell_id)));
        if isempty(sig_cells)
            sig_cells = find(cell_indices==get_cell_indices(datarun,cell_id));
            fprintf(['WARNING!!! sig_cells\n'])
        end
        
        % translate to cell indices for the whole datarun
        %sig_cells = cell_indices(sig_cells);

        
        % show stuff
        %fprintf('cell id %d (%d): %d\n',datarun.cell_ids(mm),mm,length(sig_cells))
        %figure(10);clf;hist(test_vals,500)
        %hold on;plot([1 1]*params.thresh*robust_std(test_vals),ylim,'r')
        %pause

        % create sig_stixels entry
        temp=sparse(datarun.stimulus.field_height,datarun.stimulus.field_width);
        temp(round_loc(2),round_loc(1))=1;
        for ss=1:length(sig_cells)
            sig_stixels(:,sig_cells(ss)) = reshape(temp,[],1)+sig_stixels(:,sig_cells(ss));
        end
    end
end



% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end




