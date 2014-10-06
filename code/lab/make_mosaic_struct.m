function mosaic = make_mosaic_struct(points, varargin)
% make_mosaic_struct     Compute various things for the mosaic struct
%
% usage:  mosaic = make_mosaic_struct(points, varargin)
%
% arguments:   points - Nx2 matrix of center points
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     mosaic - mosaic struct with fields
%
%      mosaic.points       - Nx2 matrix, copy of center points
%            .nnd          - Nx1 matrix, distance to each point's nearest neighbor
%            .delaunay     - Tx3 matrix, delaunay triangulation (see matlab function delaunay)
%            .nbr_dist     - NxN sparse matrix, indicates distance between delaunay neighbors
%            .delaunay_nbr - binary TxT matrix, indicates which delaunay triangles connect true neighbors
%
%
% optional params, their default values, and what they specify:
%
% verbose           false               show output
% foa               []                  figure or axes to plot in. if 0, make new figure. if empty, don't plot
% neighbor_fcn      @(nnd)(1.7*median(nnd))
%                                       cutoff for delaunay neighbors to be considered true neighbors
%
%
% 2009-02 gauthier
% 2009-09 gauthier
%
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('foa', []);
p.addParamValue('neighbor_fcn', @(nnd)(1.7*median(nnd)),@(f)isa(f,'function_handle'));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.foa);

% show output
if params.verbose
    fprintf('\nComputing mosaic properties...');
    start_time = clock; % note when it started
end


% compute NND (nearest neighbor distances)
%dists = reshape(nbr_dist(nbr_dist>0),[],1);'
dists = ipdm(points,'subset','nearest')';
nnd = reshape(dists(dists>0),[],1);


% compute delaunay triangulation
tri = delaunay(points(:,1),points(:,2));


% creat matrix of distances between delaunay neighbors

% initialize matrix
nbr_dist = sparse(size(points,1));

% go through each triangle
for tt =1:size(tri,1)
    % go through each side
    for ss=1:3
        % get point indices
        ii=tri(tt,ss);
        jj=tri(tt,mod(ss,3)+1);
        % enter these neighbor distances
        nbr_dist(ii,jj) = norm(points(ii,:)-points(jj,:));
        nbr_dist(jj,ii) = nbr_dist(ii,jj);
    end
end




% decide which triangles connect contiguous cells

% identify cutoff distance
cutoff = params.neighbor_fcn(nnd);

% note triangles with sidelengths within the cutoff
tri_nbr = false(size(tri,1));
for tt=1:size(tri,1)
    % find largest sidelength of each triangle
    largest_side = max(max(nbr_dist(tri(tt,:),tri(tt,:))));
    if largest_side <= cutoff
        tri_nbr(tt) = true;
    end
end



% plot triangulation, if desired
if ~isempty(plot_axes)
    axes(plot_axes)
    hold on
    triplot(tri,points(:,1),points(:,2),'Color',[.8 .8 .8])
    triplot(tri(tri_nbr,:),points(:,1),points(:,2),'Color',[0 0 1])
    axis ij equal
    % plot distance to nearest neighbor
    %for cc=1:length(nnd); text(points(cc,1),points(cc,2),sprintf('%0.1f',nnd(cc))); end
end


% enter values into mosaic struct
mosaic = struct;
mosaic.points = points;
mosaic.nnd = nnd;
mosaic.delaunay = uint16(tri);
mosaic.nbr_dist = nbr_dist;
mosaic.delaunay_nbr = sparse(tri_nbr);


% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end

