function returned_datarun = make_artifactual_clumping(datarun, clumped_fraction, varargin)
% make_artifactual_clumping     introduces clumping into a mosaic
%
% usage:  make_artifactual_clumping(datarun, clumped_fraction, varargin)
%
% arguments:            datarun - datarun struct with fields stimulus.stixel_height, width
%              clumped_fraction - fraction of nearest cones that get assigned the same type
%
% outputs:              returned_datarun - a new datarun is returned
%
%
% optional fields in params, their default values, and what they specify:
%
%
%
% 2009-02 gdf
%

% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addRequired('datarun', @isstruct);
p.addRequired('clumped_fraction', @isnumeric);
% specify list of optional parameters
p.addParamValue('verbose', false, @islogical);

% resolve user input and default values
p.parse(datarun, clumped_fraction, varargin{:});

verbose = p.Results.verbose;

cone_types = datarun.cones.types;
l_proportion = find(datarun.cones.types == 'L') ./ length(datarun.cones.centers(:,1));

%artifactual clumping
num_combined_cones = round(length(datarun.cones.centers(:,1)) * clumped_fraction);
temp_cone_dists = ipdm(datarun.cones.centers,'Subset', 'smallestfew', 'limit', num_combined_cones, 'result', 'struct');
flips = rand(num_combined_cones, 1);
num_changed_cones = 0;

if verbose
    figure
    hold on
end

for cone = 1:num_combined_cones
    if cone_types(temp_cone_dists.rowindex(cone)) ~= cone_types(temp_cone_dists.columnindex(cone))
        num_changed_cones = num_changed_cones + 1;
        if flips(cone) < 0.l_proportion
            cone_types(temp_cone_dists.rowindex(cone)) = 'L';
            cone_types(temp_cone_dists.columnindex(cone)) = 'L';
            
            if verbose
               plot(cone_locations(temp_cone_dists.rowindex(cone),1), cone_locations(temp_cone_dists.rowindex(cone),2), 'w.')
               plot(cone_locations(temp_cone_dists.columnindex(cone), 1), cone_locations(temp_cone_dists.columnindex(cone), 2),'w.')
               drawnow
            end
        else
            cone_types(temp_cone_dists.rowindex(cone)) = 'M';
            cone_types(temp_cone_dists.columnindex(cone)) = 'M';

            if verbose
               plot(cone_locations(temp_cone_dists.rowindex(cone),1), cone_locations(temp_cone_dists.rowindex(cone),2), 'y.')
               plot(cone_locations(temp_cone_dists.columnindex(cone), 1), cone_locations(temp_cone_dists.columnindex(cone), 2),'y.')
               drawnow
            end
        end
    end
end

returned_datarun.cones.num_changed_cones = num_changed_cones;
returned_datarun.cones.types = cone_types;

