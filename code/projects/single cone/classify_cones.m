function [cone_types,likelihood_ratios,extras] = classify_cones(cone_rgb, cone_rgb_expected, varargin)
% classify_cones     classify cones as L, M, S, or unknown
%
% usage:  [cone_types,likelihood_ratios,extras] = classify_cones(cone_rgb, cone_rgb_expected, varargin)
%
% arguments:             cone_rgb - Nx3 matrix, RGB triplet for each cone
%               cone_rgb_expected - struct with fields 'L','M','S', giving expected RGB values for each cone type
%                        varargin - struct or list of optional arguments (see below)
%
% outputs:      cone_types - Nx1 char vector, indicating the color for each cone ('L','M','S', or 'U' for unknown)
%        cone_rgb_observed - struct with fields 'L','M','S', giving observed RGB values for each cone type
%        likelihood_ratios - Nx1 matrix, likelihood ratio for each cone
%
%
% optional arguments, their default values, and what they specify:
%
%
% algorithm             'k means'
%                                   which algorithm to use to classify cones
%                                       'nearest line'
%                                       'k means'
%                                       'k means, EM'
%
%
% gauthier 2008-10
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('algorithm', 'k means, EM');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% CLASSIFY CONES


% normalize cone weights
cone_rgb_expected.L = cone_rgb_expected.L/norm(cone_rgb_expected.L);
cone_rgb_expected.M = cone_rgb_expected.M/norm(cone_rgb_expected.M);
cone_rgb_expected.S = cone_rgb_expected.S/norm(cone_rgb_expected.S);

% get number of cones
num_cones = size(cone_rgb,1);

% classify cones, and store result in cone_index
switch params.algorithm
    
    case 'nearest line'
        % for each point, find the distance to each ideal line

        % define function returning the distance from a point to a line
        point_to_line = @(pt,v1,v2)(norm(cross(v1-v2,pt-v2))/norm(v1-v2));

        field_names = {'L','M','S'};
        for ff = 1:3
            for nn = 1:num_cones
                dists(nn,ff) = point_to_line(cone_rgb(nn,:),[0 0 0],cone_rgb_expected.(field_names{ff}));
            end
        end


        % classify based on whichever line is closest
        [closest_line,junk] = find(dists' == repmat(min(dists'),3,1));

        % note the classification
        cone_index.L = find(closest_line==1);
        cone_index.M = find(closest_line==2);
        cone_index.S = find(closest_line==3);
        cone_index.U = [];
        
        % set them all to 1
        likelihood_ratios = ones(num_cones,1);

        
    case 'k means'
        % k means clustering on normalized triplets
        % initial clusters based on expected cone weights
        
        %normalize rgb values
        cone_rgb_norm = normalize_to_unit_sphere(cone_rgb);

        % run kmeans
        % assn = vector of 1s,2s,3s, giving cluster assignments
        % means = mean rgb values each cluster, one mean per row
        [assn,means] = kmeans(cone_rgb_norm,3,'start',...
            [cone_rgb_expected.L; cone_rgb_expected.M; cone_rgb_expected.S]);
        
        % identify which mean is from which cone type
        cone_types = identify_types(cone_rgb_expected,means);

        % note the classification
        cone_index.L = find(assn == cone_types.L);
        cone_index.M = find(assn == cone_types.M);
        cone_index.S = find(assn == cone_types.S);
        cone_index.U = [];
        
        % note likelihood ratios
        likelihood_ratios = ones(num_cones,1);
        
        
    case 'k means, EM'
        % k means clustering on normalized triplets
        % initial clusters based on expected cone weights
        
        %
        % parameter for rejecting points (in likelihood ratios):
        %

        threshold = 0;

        %
        % step one: assign points to clusters
        %

        % fix the seed value
        rand('state',11111)

        [typeKmeans, means] = kmeans(cone_rgb, 3, 'EmptyAction', 'drop', 'Start',  ...
            [cone_rgb_expected.L; cone_rgb_expected.M; cone_rgb_expected.S],...
            'Distance', 'cosine');

        %
        % step two: reject points where assignments are uncertain
        %

        %fit gaussians and compute likelihood ratios
        nDims = 3;
        nClusters = 3;

        %normalize points to the unit sphere
        cn = zeros(length(cone_rgb),3); meansN = zeros(size(means));
        for ii=1:length(cone_rgb)
            cn(ii,:) = cone_rgb(ii,:)./norm(cone_rgb(ii,:));

        end

        %also normalize centroids
        for ii=1:nClusters
            meansN(ii,:) = means(ii,:)./norm(means(ii,:));
        end

        %parameters needed by by vision (these are the defaults)
        likelihoodDelta = 1e-6;
        minIterations   = 5;
        maxIterations   = 200;

        %set up em: we just use this code to compute the likelihood of the
        %gaussians that we manually specify based on the output from kmeans
        em = edu.ucsc.neurobiology.vision.math.ExpectationMaximization1(nDims, nClusters, length(cn));
        em.reset(cn',length(cn),length(cn));

        %don't allow em to fit the mean
        em.setState(false, true, true);

        for ii = 1:nClusters
            %use the variance of the cluster from kmeans
            v = var(cn(typeKmeans==ii,:));

            %use means from kmeans
            em.addGaussian(meansN(ii,:),v);
        end

        %fit the clusters from kmeans using em
        em.fit(likelihoodDelta, minIterations, maxIterations);

        %compute the likelihood for all points
        constant = (2*pi)^(nDims/2);
        likelihood = zeros(length(cn), nClusters);
        typeEM    = zeros(length(cn), 1);
        lRatio     = zeros(length(cn), 1);

        for point=1:length(cn)
            %compute the likelihood of each point given the optimized gaussian fits
            for cluster=0:nClusters-1
                likelihood(point,cluster+1) = javaMethod('_lpXJ', em, cn(point,:),em.getMeans(cluster), em.getCovariances(cluster), constant);
            end

            li = sort(likelihood(point,:));
            typeEM(point) = find(likelihood(point,:) == li(end), 1, 'last' );
            lRatio(point) = li(end) - li(end-1);
        end
        
        
        % output variables
        %
        % typeKmeans    class assignments from kmeans
        % typeEM        class assignments with points rejected at threshold
        % lRatio        likelihood ratios for all points

        
        
        
        % identify points that changed class, label them as unknown
        
        % initialize to EM classification
        type = typeEM;
        % in each cell...
        for point=1:length(cn)
            % if the assignment changed, reject the point
            if(typeEM(point) ~= typeKmeans(point))
                type(point) = 4;
            end
        end
        
        % get mean of each class
        for mm=1:4
            means(mm,:) = mean(cn(type==mm,:));
        end
        
        % identify which mean is from which cone type
        cone_types = identify_types(cone_rgb_expected,means);

        % note the classification
        cone_index.L = find(type == cone_types.L);
        cone_index.M = find(type == cone_types.M);
        cone_index.S = find(type == cone_types.S);
        cone_index.U = find(type == cone_types.U);
        
        % note likelihood ratios
        likelihood_ratios = lRatio;
        
        % enter other classification in to extras
        extras.types_em = cone_type_nums_to_chars(typeEM);
        extras.types_kmeans = cone_type_nums_to_chars(typeKmeans);
        
        
        
    otherwise
        error('Algorithm %s not recognized.',params.algorithm)
        
end



% set cone_types variable

% initialize
cone_types = char(1,num_cones);

% make assignments based on cone_index
cone_types(cone_index.L) = 'L';
cone_types(cone_index.M) = 'M';
cone_types(cone_index.S) = 'S';
cone_types(cone_index.U) = 'U';


% compute mean RGB values for each cone type
cone_rgb_observed.L = mean(cone_rgb(cone_types == 'L',:));
cone_rgb_observed.M = mean(cone_rgb(cone_types == 'M',:));
cone_rgb_observed.S = mean(cone_rgb(cone_types == 'S',:));



% set up struct of extras
extras.cone_rgb_observed = cone_rgb_observed;







function cone_types = identify_types(cone_rgb_expected,means)
% match the means of cone clusters to their identiy (L,M,or S) based on expected cone weights
% means - mean rgb values each cluster, one mean per row
%
% returns a struct with fields L, M, S, U, indicating which row of means corresponds to which cone type
%       e.g. if cone_types.L == 2, then the second row of means corresponds to the L cones
%

% normalize to the unit sphere
means = normalize_to_unit_sphere(means);


% get list of possible re-assignments
ra = perms(1:3);

% compute distances for each reassignment
for rr=1:size(ra,1)
    % make new means variable with reassignments
    means_temp = means(ra(rr,:),:);
    
    % get "error" of this assignment as the sum of squared distances 
    err(rr) = norm([cone_rgb_expected.L; cone_rgb_expected.M; cone_rgb_expected.S] - means_temp);
    
end

% find which possible re-assignment had the smallest error
[junk,rr_best] = min(err);

% make assignments based on this permutation
reassigner = ra(rr_best,:);
cone_types.L = reassigner(1);
cone_types.M = reassigner(2);
cone_types.S = reassigner(3);
% 4 is always unknown
cone_types.U = 4;




function chars = cone_type_nums_to_chars(nums)

chars = char;
chars(nums == 1,1) = 'L';
chars(nums == 2,1) = 'M';
chars(nums == 3,1) = 'S';
chars(nums == 4,1) = 'U';

