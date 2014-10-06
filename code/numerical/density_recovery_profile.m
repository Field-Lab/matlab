function [drp,bin_centers,extras] = density_recovery_profile(centers,num_bins,bin_width,varargin)
% density_recovery_profile     Compute density recovery profile and related tidbits
%
%    Everything is based on this paper:
%
%   Rodieck, RW "The density recovery profile: A method for the analysis of points in the plane applicable to retinal studies"
%       Visual Neuroscience (1991) 6, 95-111
%
%
% usage:  [drp,bin_centers,extras] = density_recovery_profile(centers,num_bins,bin_width,varargin)
%
% arguments:     centers - Nx2 matrix, coordinates of center points
%               num_bins - how many bins
%              bin_width - 
%               varargin - struct or list of optional parameters (see below)
%
% outputs:         drp - estimated density in each annulus
%
%
% optional params, their default values, and what they specify:
%
% reference_centers     []      Mx2 matrix, coordinates of points to use as reference for the DRP
% foa                   []      figure or axes where to plot a bar graph of the results.  if empty, don't plot
%
%
% 2008-10 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('reference_centers', []);
p.addParamValue('foa', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


if isempty(params.reference_centers)
    % get all pair-wise distances

    switch 1
        case 1
            % use the points themselves
            dists = pdist(centers);

            % normalization for number of points
            % for an auto-correlation, this is simply the number of points
            num_norm = size(centers,1);

            % also, should divide by 0.5.  under rodieck's method, each distance should be counted twice
            % num_norm is applied as a divisive normalization below
            % thus multiplying by 0.5 doubles the count in each bin, effectively
            % counting every distance twice
            num_norm = 0.5*num_norm;
            
            
        case 2
            dists = ipdm(centers);%,'Subset','Max','limit',num_bins*bin_width);
            dists = squareform(dists);
            num_norm = size(centers,1);

    end
else
    
    % otherwise, get distances from the reference centers to the other centers

    switch 2
        case 1
            % combine all centers, with reference centers first
            all_centers = [params.reference_centers; centers];

            % get all distances
            all_dists = pdist(all_centers);

            % make full distance matrix
            dist_mat = squareform(all_dists);

            % get upper right section, which correspondes to distances between reference and non-reference points
            dists = reshape(dist_mat(1:size(params.reference_centers,1),end-size(centers,1)+1:end),[],1);

        case 2
            dists = ipdm(params.reference_centers,centers,'Subset','Max','limit',num_bins*bin_width);
            dists = reshape(dists,[],1);

    end


    % note number of points
    % for an auto-correlation, this is geometric mean of the number of reference points and the number of other points
    %num_norm = sqrt(size(params.reference_centers,1)*size(centers,1));
    num_norm = size(params.reference_centers,1);
end

% get bin centers, adding an extra bin at the end to catch points which are farther away.
% this extra bin will be discarded later
bin_centers = bin_width/2 + (0:1:num_bins)*bin_width;

% bin distances
counts = hist(dists(dists<Inf),bin_centers);

% discard the last bin
counts = counts(1:end-1);

% get area of each annulus
areas = pi*bin_width^2*(2*(1:num_bins)-1);

% divide by the area and the normalization for the number of points (defined above)
drp = counts./(areas*num_norm);


% for output, discard the last bin
bin_centers = bin_centers(1:end-1);


% get additional info, if desired
if nargout>2
    
    % mean density
    % average over the larger half of the bins
    extras.density = mean(drp(round(end/2):end));
    
    % effective radius
    % identify the bins with density >= mean density
    temp = find(drp >= extras.density);
    % if the first bin has this density, then the effective radius can not be estimated
    if temp(1) == 1
        extras.eff_rad = 0;
    else
        % otherwise, get the index of the furthest bin with density < mean density
        furthest_bin = temp(1) - 1;
        
        % compute the total missing "volume" of the bins
        %eff_vol = sum([extras.density - drp(1:furthest_bin)].*areas(1:furthest_bin))/furthest_bin;
        
        % expected counts
        expected_counts = extras.density .* (areas*num_norm);
        
        % compute the total missing "volume" of the bins
        eff_vol = sum(expected_counts(1:furthest_bin) - counts(1:furthest_bin))/num_norm;
        
        % compute the effective radius
        extras.eff_rad = sqrt(eff_vol/(pi*extras.density));
    end
    
    
    
    %     % effective radius
    %     % identify the furthest bin with density < mean density
    %     temp = find(drp < extras.density);
    %     furthest_bin = temp(end);
    %     % compute the total missing "volume" of the bins
    %     eff_vol = sum(drp(1:furthest_bin));
    %     % compute the effective radius
    %     extras.eff_rad = sqrt(eff_vol/(pi*extras.density));

end



% plot, if desired
if ~isempty(params.foa)
    plot_axes = set_up_fig_or_axes(params.foa);
    axes(plot_axes)
    bar(bin_centers,drp)
end
