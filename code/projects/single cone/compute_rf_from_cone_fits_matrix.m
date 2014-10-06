function M = compute_rf_from_cone_fits_matrix(cone_centers,cone_rgb,varargin)
% compute_rf_from_cone_fits_matrix     compute matrix so that
%
%   rf_fits = M * cone_weights
%
%  where
%
%       cone_weights - NxR matrix, N = # cones, R = # RFs
%            rf_fits - vector such that reshape(rf_fits,y_size,x_size,3,num_rfs)
%                       contains RF fits for the given cone weights
%
% usage:  M = compute_rf_from_cone_fits_matrix(cone_centers,cone_rgb, varargin)
%
% arguments:     cone_ctrs - Nx2 matrix
%                 cone_rgb - Nx3 matrix
%                 varargin - struct or list of optional parameters (see below)
%
% outputs:     M - matrix, see above
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
% 2009-04  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('x_size', 10);
p.addParamValue('y_size', 10);
p.addParamValue('num_rfs', 3);
p.addParamValue('center_radius', 0.75);
p.addParamValue('fig', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION



M = sparse(0,0);

for cc = 1:size(cone_centers,1)
    
    % make bw kernel
    bw_kernel = make_gaussian(struct('center',cone_centers(cc,1:2),...
        'y_size',params.y_size,'x_size',params.x_size,...
        'center_radius',params.center_radius,'normalize','sum',...
        'effective_radius',round(3*params.center_radius)));

    % make RGB kernel
    rgb_kernel = zeros(params.y_size,params.x_size,size(cone_rgb,2));
    for dd = 1:size(cone_rgb,2)
        % weight cone kernel for this color
        rgb_kernel(:,:,dd) = bw_kernel*cone_rgb(cc,dd);
    end
    rgb_kernel = sparse(reshape(rgb_kernel,[],1));


    % add to accumulating matrix
    temp = cell(params.num_rfs,1);
    for rr=1:params.num_rfs
        temp{rr}=rgb_kernel;
    end
    M = [M blkdiag(temp{:})];

end





if ~isempty(params.fig)

    figure(params.fig);clf;
    plot_axes = subplot_axes(params.fig,[0 0 1 1],0.1,0.1,params.num_rfs,1);

    %cw = [1 0 0 0  0 1 0 0  0 0 1 1]';
    %cw = [1 0 0  1 0 0  0 1 1  0 0 1]';
    cw = [1 0 0;...
          0 1 0;...
          0 1 1;...
          0 0 0];
    
    cw = reshape(cw',[],1);

    rfs = reshape(M*cw,params.y_size,params.x_size,3,params.num_rfs);
    
    reshape(cw,params.num_rfs,size(cone_centers,1))'

    for rr=1:size(rfs,4)
        axes(plot_axes{rr})
        imagesc(norm_image(rfs(:,:,:,rr)))
        axis image;hold on
        for cc=1:size(cone_centers,1)
            plot(cone_centers(cc,1),cone_centers(cc,2),'.')
            text(cone_centers(cc,1),cone_centers(cc,2)-0.5,sprintf('cone %d',cc))
        end
    end

end

