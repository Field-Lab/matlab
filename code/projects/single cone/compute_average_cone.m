function avg_cone = compute_average_cone(cone_rfs, varargin)
% MY_FUNCTION     This template function does nothing.  Ha!
%
% usage:  result = my_function(arg1, <params>)
%
% arguments:     arg1 - argument 1
%              params - struct or list of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional params, their default values, and what they specify:
%
% verbose           false               show output
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure. if empty, don't plot
% foo               'bar'               how to activate the rings
%                                           'bar' - activate on site
%                                           'bore' - activate remotely
%
%
% author date
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional arguments
p.addParamValue('verbose', false);
p.addParamValue('fig_or_axes', []);
p.addParamValue('foo','bar', @(x)any(strcmpi(x,{'bar','bore'})));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% initialize variable
avg_cone = zeros(3,3,3);

% get each observed cone rf
for nn =1:length(cone_rfs)
    big_cone_rf(:,:,1) = full(cone_rfs{nn}.cone_map_sum_r);
    big_cone_rf(:,:,2) = full(cone_rfs{nn}.cone_map_sum_g);
    big_cone_rf(:,:,3) = full(cone_rfs{nn}.cone_map_sum_b);

    % find non-zero points
    [ii,jj] = find(sum(big_cone_rf,3) ~= 0);

    % get them
    cone_neighborhood = big_cone_rf(unique(ii),unique(jj),:);

    % if they are 3x3, add them
    if all(size(cone_neighborhood) == [3 3 3])

        % enter values from the cone region
        avg_cone(:,:,:) = avg_cone(:,:,:) + cone_neighborhood;

    end

    if ~isempty(plot_axes)
        axes(plot_axes)
        cla
        imagesc(norm_image(avg_cone))
        axis image
        drawnow
    end
end




