function cone_roi = identify_cone_mosaic_roi(centers, cone_density, cone_eff_rad)
% identify_cone_mosaic_roi     identify a ROI in the cones based on apparently contiguous regions
%
% usage:  cone_roi = identify_cone_mosaic_roi(centers, cone_density, cone_eff_rad)
%
% arguments:  centers - datarun struct with field specifying X
%        cone_density - density of cones 
%        cone_eff_rad - effective radius of the cone mosaic
%
% outputs:    cone_roi - binary vector
%
%
% 2008-10 gauthier
%




% BODY OF THE FUNCTION

switch 2
    case 1 % neighbors
        % apply delaunay triangulation
        tri = delaunay(centers(:,1),centers(:,2));

        % get area and longest sidelength of each triangle
        tri_area = zeros(size(tri,1),1);
        tri_max_side = zeros(size(tri,1),1);
        for tt=1:size(tri,1)
            vx = centers(tri(tt,:),1);
            vy = centers(tri(tt,:),2);
            tri_area(tt) = polyarea(vx,vy);
            tri_max_side(tt) = max(pdist([vx vy]));
        end

        % identify triangle sconnecting "adjacent" cones.
        good_tri = tri_area<Inf*(1/density) & tri_max_side < 3*eff_rad;

        % identify points which border at least one such triangle
        nbr_cones = unique(tri(good_tri,:));

        cone_roi = false(size(centers,1),1);
        cone_roi(nbr_cones) = true;
        
        
    case 2
        
        % make a filter to count number of cones within a certain radius
        
        % inclusion radius
        radius = 5*cone_eff_rad;
        
        % make disk of this radius
        filt_size = round([1 1]*(2*radius+1));
        filt = double(distance_from_point(filt_size,[1 1]*(radius + 1)) <= radius);
        
        
        % make a matrix with each cone located at its rounded coordinates, 
        % and label each with a different number
        cone_locs = cone_rf_reconstructed(round([max(centers(:,2)) max(centers(:,1))]+radius),1:size(centers,1),centers);
        
        
        % use the filter to get local density at every point
        local_density = imfilter(double(cone_locs>0),filt);
        
        
        % set density threshold
        thresh = 1.1*sum(sum(filt))*cone_density;
        %thresh = 0.4*sum(sum(filt))*(0.5/(cone_eff_rad^2));
        
        
        % identify cones centered on a location where the local density exceeds threshold
        cone_indices = reshape((local_density>thresh).*cone_locs,[],1);
        cone_indices = setdiff(cone_indices,0);

        % make the list of cone indices into a binary vector
        cone_roi = false(size(centers,1),1);
        cone_roi(cone_indices) = true;
        
        % plot it
        %figure(9);clf;image(1*(local_density>thresh)+32);axis image;colormap('gray')
        %hold on; plot(centers(:,1),centers(:,2),'.k',centers(cone_roi,1),centers(cone_roi,2),'.b')
end

