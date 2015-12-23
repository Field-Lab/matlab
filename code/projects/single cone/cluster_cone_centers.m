function [clusters,cluster_centers,cluster_sources,max_dists] = cluster_cone_centers(datarun,points,cone_sources,varargin)
% cluster_cone_centers     Combine local maxima found in each RF to identify which are from the same cone
%
% usage:  datarun = my_function(datarun, arg1, <params>)
%
% arguments:  datarun - datarun struct with field specifying X
%        points - argument 1
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in field X
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           true                show output
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure. if empty, don't plot
% radius_1          1.2
%
%
%
% date author
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', true);
p.addParamValue('radius_1', 1.2);
p.addParamValue('radius_2', 1.6);
p.addParamValue('fig_or_axes', []);


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);


% cluster to find connected points
T=cluster(linkage(pdist(points)),'cutoff',params.radius_1,'criterion','distance');

% get number of clusters
num_clust = length(unique(T));

% show output
if params.verbose
    fprintf('Found %d clusters from %d local maxima\n',num_clust,size(points,1));
end

% find largest distance in each
[max_dist,tt] = max(md(points,T));

% if it's too big, subcluster
while max_dist > params.radius_2
    
    fprintf('cluster %d has max distance %0.3f\n',tt,max_dist)

    % get points from largest cluster
    these_points = points(T==tt,:);
    
    % perform linkage
    Z = linkage(pdist(these_points));
    
    % cluster at the level of the highest break
    U = cluster(Z,'cutoff',Z(end,end)-.000001,'criterion','distance');
    
    % insert into list of clusters
    % U==1  ->  tt
    % U==2  ->  num_clusts + 1
    % U==3  ->  num_clusts + 2
    % ...
    V=U;
    V(U==1)=tt;
    V(U>1)=U(U>1)+length(unique(T))-1;
    T(T==tt)=V;
    
    % identify new largest distance
    [max_dist,tt] = max(md(points,T));
    
end


% get number of clusters
num_clust = length(unique(T));

% find center and sources for each cluster
cluster_centers = zeros(num_clust,2);
cluster_sources = cell(num_clust,1);
for tt=1:num_clust
    % take mean of individual centers
    cluster_centers(tt,:) = mean(points(T==tt,:),1);
    % list all sources
    cluster_sources{tt} = cone_sources(T==tt);
end

% get maximum distance of points in each cluster
max_dists = md(points,T);

% change to output name
clusters = T;



% find maximum distance of each cluster
function [max_dists] = md(points,T)


% get number of clusters
num_clust = length(unique(T));

% initialize max distances to 0
max_dists=zeros(num_clust,1);

% get actual max distance of each cluster
for tt=1:num_clust
    dists = pdist(points(T==tt,:));
    if isempty(dists)
        max_dists(tt) = 0;
    else
        max_dists(tt) = max(dists);
    end
end
