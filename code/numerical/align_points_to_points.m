function remap = align_points_to_points(orig_points,ref_points,varargin)
% align_points_to_points     align a set of points and to a different set of reference points
%                               a transformation may also be applied to the points
%
% usage:  point_tform = align_points_to_points(orig_points,ref_points,varargin)
%
% arguments:     orig_points - Nx2 matrix giving coordinates for the points to be mapped
%                 ref_points - Nx2 matrix giving coordinates for reference points
%
% outputs:     remap - NxM matrix mapping the original points to the reference points
%
%
%
% optional params, their default values, and what they specify:
%
% tform             []              optional transformation to apply
%                                       must be in the format of a "spatial transformation structure"
%                                       see "maketform" and related functions
% fig_or_axes       []             	where to plot the points before and after
%                                      	if 0, make new figure. if empty, don't plot
%
%
% 2009-02  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('tform', []);
p.addParamValue('verbose', 0);
p.addParamValue('method', 'nearest', @(x)any(strcmpi(x,{'nearest'})));
p.addParamValue('fig_or_axes', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% show output
if params.verbose
    fprintf('\nComputing something important...');
    start_time = clock; % note when it started
end

% apply optional transformation
if ~isempty(params.tform)
    orig_points = tformfwd(orig_points,params.tform);
end
    

% align points
switch params.method
    case 'nearest'
        
        % for each original point, get get nearest reference point
        d=ipdm(orig_points,ref_points,'Subset','nearest');
        % make into remap matrix
        remap = full(d)>0;
        
    otherwise
        error('Method not recognized: %s',params.method)
end



if ~isempty(plot_axes)
    axes(plot_axes); hold on
    % plot the points
    plot(orig_points(:,1),orig_points(:,2),'r.','MarkerSize',20)
    plot(ref_points(:,1),ref_points(:,2),'k.','MarkerSize',20)
    % plot lines connecting each to its nearest neighbor
    for pp = 1:size(orig_points,1)
        rr = find(remap(pp,:));
        plot([orig_points(pp,1) ref_points(rr,1)],[orig_points(pp,2) ref_points(rr,2)],'b')
    end
end


% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end

