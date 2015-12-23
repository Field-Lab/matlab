function new_selection = enforce_cone_contiguity(datarun, weights, selection, scale)
% enforce_cone_contiguity     given a set of cones and weights, identify which cones are contiguous with the strongest cone
%
% usage:  new_selection = enforce_cone_contiguity(datarun, weights, selection, scale)
%
% arguments:     datarun - datarun struct
%                weights - double NxM vector of cone weights
%              selection - logical NxM vector of which cones to choose from
%                  scale - scale factor for determining whether cones are contiguous
%                             cutoff = scale * median cone nearest neighbor spacing
%
% outputs:     new_selection - logical N-vector of which cones are contiguous
%
%
% 2009-10 gauthier, based on procedure developed by GDF
%   2009-10 gdf, introduced for loop to enable function to handle many RGCs
%   2010-4 gdf, fixed bug that occured when connectivity was only to 1 cone
%



if 1
    [num_cones, num_rgcs] = size(selection);

    % initialize output
    new_selection = false(size(selection));

    for cc = 1:num_rgcs

        % get the selection vector for one rgc
        temp_selection = selection(:,cc);
        temp_weights = weights(:,cc);

        % get points
        points = datarun.cones.centers(temp_selection,:);


        if size(points,1) == 1
            new_selection(:,cc) = temp_selection;

        else

            % get index of strongest cone
            [junk,strongest_index] = max(temp_weights(temp_selection));

            % cluster points
            T = cluster(linkage(pdist(points)),'cutoff',scale * median(datarun.cones.mosaic.nnd),'criterion','distance');

            % identify cluster of strongest index
            strongest_cluster = T(strongest_index);

            keep_indices = find(T == strongest_cluster);

            % temp_selection is points contiguous with this point
            new_selection(temp_selection,cc) = (T == strongest_cluster);

        end
    end
end




% previous code
if 0


    % created: 2009-07-13, GDF

    % body of function


    % get number of RGCs and number of cones
    [num_cones, num_RGCs] = size(weights);

    % initialize new_selection
    new_selection = false(num_cones, num_RGCs);

    connectivity = weights .* selection;

    % get median cone nnd
    %nnds = ipdm(datarun.cones.centers, 'Subset', 'NearestNeighbor');
    %nnd_indices = find(nnds >0);
    %nnds = nnds(nnd_indices);
    %median_nnd = median(nnds);

    radius = scale * median(datarun.cones.mosaic.nnd);


    for RGC = 1:num_RGCs


        % for an RGC, get the non-zero weighted (selected) cones
        temp_cone_indices = find(connectivity(:, RGC) ~= 0);

        if ~isempty(temp_cone_indices)


            [sorted_weights, sorted_indices] = sort(connectivity(temp_cone_indices, RGC), 'descend');

            % the cone indices are sorted by cone weights
            sorted_cone_indices = temp_cone_indices(sorted_indices);

            % seed the center cone locations with the location of strongest cone
            center_locations = datarun.cones.centers(sorted_cone_indices(1),:);
            tracker_indices = sorted_cone_indices(1);

            for cone = 1:num_cones

                % indices to cones not yet included in the center
                [remaining_cone_indices, remaining_indices] = setdiff(sorted_cone_indices, tracker_indices);

                % find distances between "center cones" and remaining cones
                temp_distances = ipdm(center_locations, datarun.cones.centers(sorted_cone_indices(remaining_indices),:));


                % find the strongest weighted cone that is contiguous
                temp_index = find(temp_distances < radius, 1);

                % temp_c is the index to the strongest cone that satisfied the
                % contiguity criterion
                [temp_r, temp_c] = ind2sub(size(temp_distances), temp_index);

                % update indices to contiguous cones
                tracker_indices = [tracker_indices, remaining_cone_indices(temp_c)];

                % add new cone location to center_locations
                center_locations = [center_locations; datarun.cones.centers(remaining_cone_indices(temp_c),:)];

                % if no more cones or all cones satisfy contiguity, then stop
                if isempty(temp_index) || length(tracker_indices) == length(temp_cone_indices)
                    break
                end

            end

            % generate a selection matrix for contiguous cones
            contig_selection = false(num_cones,1);
            new_selection(tracker_indices,RGC) = true(length(tracker_indices),1);
            
        end

    end

end

