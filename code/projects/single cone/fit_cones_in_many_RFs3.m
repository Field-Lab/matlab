function [fit_cone_info,rf_fits] = fit_cones_in_many_RFs3(rfs,cone_info,varargin)
% FIT_CONES     fit ideal cone shapes to cones in several RFs
%
% #3 removed input strengths
%
% usage:  [fit_cone_info,rf_fit] = fit_cones(rfs,cone_locations,params)
%
% arguments:       rfs - YxXxCxM matrix of intensity map of each RF, C = # colors, M = # RFs
%            cone_info - Nx3 matrix, initial guesses of cone locations, N = # cones
%                           column 1: x coord
%                           column 2: y coord
%                           column 3: cone identity
%
% output:   fit_cone_info - Nx(M+3) matrix, fits to cone locations, some format as cone_info
%                 rf_fits - YxXxC,M matrix of fits to each RF
%
% optional fields in params, their default values, and what they specify:
%
% figure            []          which figure to plot in while being fit.  if 0, make new.  if empty, don't plot
% fig_final         []          where to plot the final fit
%
% cone_kernel       struct('type','dog')
%                               what kind of shape to fit
%
%                               struct of parameters to generate the cone kernel which are regressed against the cells
%
%                                   type        how to generate cone kernel
%                                                   'dog'     - difference of gaussian
%
%
%                                   only required field is 'type'
%                                   all other fields are optional, with different default values depending on the type of kernel
%
%                                   for type = 'dog', normalization is set to 'sum', and
%                                       these are the optional fields (see 'make_gaussian'):
%
%                                           center_radius       0.75    center gaussian radius
%                                           surround_radius     2.5     surround gaussian radius
%                                           surround_scale      0       surround gaussian scale (NOT relative)
%                                           effective_radius    2       region in which to compute the gaussian
%
%
% cone_rgb          struct with fields
%                       L: [0.4044 0.8854 0.2292]
%                       M: [0.1483 0.9161 0.3726]
%                       S: [0.0174 0.0792 0.9967]
%                                   expected RGB senstivities for for the cones of each type
%
% optim             {'TolX',.1}
%                               arguments to optimset that control optimization of fit.
%                                   see documentation of fminsearch for options.
%
% gauthier, 2008-09
%




% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('cone_kernel', struct('type','dog'));
p.addParamValue('figure', []);
p.addParamValue('fig_final', []);
p.addParamValue('optim',{'TolX',0.1,'Display','off'});
%p.addParamValue('algorithm', 'fminsearch',@(x)any(strcmpi(x,{'fminsearch','lsqcurvefit'})));
p.addParamValue('algorithm', 'lsqcurvefit',@(x)any(strcmpi(x,{'fminsearch','lsqcurvefit','lsqnonlin','fmincon','com'})));
p.addParamValue('cone_rgb', struct(...
    'L',[0.4044 0.8854 0.2292],...
    'M',[0.1483 0.9161 0.3726],...
    'S',[0.0174 0.0792 0.9967]));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% get default kernel parameters
clear defaults
if ~isempty(params.cone_kernel)
    switch params.cone_kernel.type

        case 'dog'
            defaults.type = 'dog';
            defaults.center_radius = 0.75;
            defaults.surround_radius = 2.5;
            defaults.surround_scale = 0;
            defaults.effective_radius = 2;

        otherwise
            error('Cone kernel type ''%s'' not valid',params.cone_kernel.type)
    end

    % combine user and default parameters
    params.cone_kernel = default_params( defaults, params.cone_kernel);
end

% add cone_rgb to cone_kernel params
params.cone_kernel.cone_rgb = params.cone_rgb;



% set up figure
if ~isempty(params.figure)
    if params.figure > 0
        figure(params.figure)
    else
        params.figure = figure;
    end
    figure(params.figure);clf
    plot_axes = subplot_axes(params.figure,[0 0 1 1],0.1,0.1,size(rfs,4),1);
end


% start with approximate cone colors
[cone_weights,cone_colors] = find_best_weights(rfs,cone_info,params.cone_kernel);
cone_info(:,end) = cone_colors;


% generate upper and lower bounds
min_x = 0.5;
min_y = 0.5;
max_x = size(rfs,2)+0.5;
max_y = size(rfs,1)+0.5;
lb = repmat([min_x min_y 0],size(cone_info,1),1);
ub = repmat([max_x max_y 1],size(cone_info,1),1);


% optimize fit, save cone info for best fit
switch params.algorithm
    case 'fminsearch'
        obj_fcn = @(x)get_fit_err(rfs,x,params.cone_kernel);
        plot_fcn = @(xx,oo,ss)plot_current_fits(xx,oo,ss,rfs,params.cone_kernel,plot_axes);

        fit_cone_info = fminsearch(obj_fcn,...
            cone_info,... initial parameters
            optimset(params.optim{:},...  % optional user parameters
            'OutputFcn',plot_fcn,... plotting function
            'Display','iter')); %

    case 'lsqcurvefit'
        % make function that reads in the parameters AND the x data
        fit_fcn = @(cone_info,xx)get_rf_fits(rfs,cone_info,params.cone_kernel);
        plot_fcn = @(xx,oo,ss)plot_current_fits(xx,oo,ss,rfs,params.cone_kernel,plot_axes);

        fit_cone_info = lsqcurvefit(...
            fit_fcn,...
            cone_info,... initial parameters
            1,... x data (not relevant here)
            rfs,... y data
            lb,ub,... upper, lower bounds
            optimset(params.optim{:},... optional user parameters
            'OutputFcn',plot_fcn)); % plotting function

    case 'lsqnonlin'
        obj_fcn = @(x)get_fit_err(rfs,x,params.cone_kernel);
        plot_fcn = @(xx,oo,ss)plot_current_fits(xx,oo,ss,rfs,params.cone_kernel,plot_axes);

        fit_cone_info = lsqnonlin(obj_fcn,...
            cone_info,... initial parameters
            lb,ub,... upper, lower bounds
            optimset(params.optim{:},...  % optional user parameters
            'OutputFcn',plot_fcn)); % plotting function

    case 'fmincon'
        obj_fcn = @(x)get_fit_err(rfs,x,params.cone_kernel);
        plot_fcn = @(xx,oo,ss)plot_current_fits(xx,oo,ss,rfs,params.cone_kernel,plot_axes);

        fit_cone_info = fmincon(obj_fcn,...
            cone_info,... initial parameters
            [],[],[],[],lb,ub,[],... upper, lower bounds
            optimset(params.optim{:},...  optional user parameters
            'OutputFcn',plot_fcn,... plotting function
            'RelLineSrchBnd',0.5,'RelLineSrchBndDuration',Inf,'algorithm','trust-region-reflective')); %

    otherwise
        error('fitting algorithm ''%s'' not recognized.',params.algorithm)
end

% return the fit with the best parameters
rf_fits = get_rf_fits(rfs,fit_cone_info,params.cone_kernel);

% plot final fit, if desired
if ~isempty(params.fig_final)

end




function err = get_fit_err(rfs,cone_info,kernel_params)
% take in cone locations and return error rate

% get RF fit based on current cone info
rf_fits = get_rf_fits(rfs,cone_info,kernel_params);

% compute error
err = (mean(reshape((rf_fits - rfs).^2,[],1)));
%err = reshape((rf_fits - rfs).^2,[],1);

% add in penalty functions
if 0

    % cone distance penalty
    dist_penalty = 0;
    dist_pen_fcn = @(x)1000*x/exp(0.7*x^2);
    cone_centers = cone_info(:,end-2:end-1);
    for cc = 1:size(cone_centers,1)
        for dd=1:cc-1
            dist = norm(cone_centers(cc,:)-cone_centers(dd,:));
            dist_penalty = dist_penalty + dist_pen_fcn(dist);
        end
    end

    err = err + dist_penalty;
end



% return a fit
function rf_fits = get_rf_fits(rfs,cone_info,kernel_params)

% more convenient variables
cone_centers = cone_info(:,1:2);
cone_colors = cone_info(:,3);

num_cones = size(cone_info,1);
cone_fit_size = [size(rfs,1) size(rfs,2)];


% get cone RGB
the_rgb = zeros(num_cones,3);
for cc = 1:num_cones
    the_rgb(cc,:) = cone_rgb(cone_colors(cc),kernel_params.cone_rgb);
end

% get matrix
M = compute_rf_from_cone_fits_matrix(cone_centers,the_rgb,'x_size',size(rfs,2),'y_size',size(rfs,1),...
    'center_radius',kernel_params.center_radius,'num_rfs',size(rfs,4));

switch 3
    case 1  % compute weights for each cone
        best_weights = M\reshape(rfs,[],1);
    case 2  % make all 1s
        best_weights = 1*ones(size(M,2),1);
    case 3  % compute weights for each cone, preserve only sign
        best_weights = 0.1*sign(M\reshape(rfs,[],1));
end

% compute RF fits in same format as RFs
rf_fits = reshape(M*best_weights,size(rfs));



% plot RFs and fits
if 1
    figure(13);clf
    plot_axes = subplot_axes(13,[0 0 1 1],0.1,0.1,size(rfs,4),1);

    % plot them
    pa = plot_axes;
    for rr = 1:size(rfs,4)

        % plot RFs
        cla(pa{rr})
        %imagesc(norm_image([rfs(:,:,:,rr);rf_fits(:,:,:,rr);rf_fits_p1(:,:,:,rr);rf_fits_p2(:,:,:,rr)]),'Parent',pa{rr});
        imagesc(norm_image([rfs(:,:,:,rr);50*rf_fits(:,:,:,rr)]),'Parent',pa{rr});axis(pa{rr},'image')
        % plot cones
        hold(pa{rr},'on')
        % over RFs
        plot(pa{rr},cone_info(:,end-2),cone_info(:,end-1),'.k')
        % over fits
        plot(pa{rr},cone_info(:,end-2),cone_info(:,end-1)+size(rfs,1),'.k')
    end
    drawnow
end




function rgb = cone_rgb(alpha,cone_rgb)
% return RGB triplet, given alpha

% constrain to be between 0 and 1
if alpha < 0; alpha = 0; end
if alpha > 1; alpha = 1; end

% take weighted sum of L/M and S
rgb = (1-alpha)*(cone_rgb.L+cone_rgb.M)/2 + (alpha)*cone_rgb.S;







function stop = plot_current_fits(x,optimValues,state,rfs,kernel_params,plot_axes)

%fprintf('plot fcn called (iteration %d, state %s, resnorm %0.1f)\n',optimValues.iteration,state,optimValues.resnorm)
fprintf('plot fcn called (iteration %d, state %s)\n',optimValues.iteration,state)

stop = false;

if mod(optimValues.iteration,1)==0
    % get current parmeters
    cone_info = x;

    % generate fits
    rf_fits = get_rf_fits(rfs,cone_info,kernel_params);

    % plot them
    pa = plot_axes;
    for rr = 1:size(rfs,4)

        % plot RFs
        cla(pa{rr})
        imagesc(norm_image([rfs(:,:,:,rr);rf_fits(:,:,:,rr)]),'Parent',pa{rr});axis(pa{rr},'image')
        % plot cones
        hold(pa{rr},'on')
        % over RFs
        plot(pa{rr},cone_info(:,end-2),cone_info(:,end-1),'.k')
        % over fits
        plot(pa{rr},cone_info(:,end-2),cone_info(:,end-1)+size(rfs,1),'.k')
    end
    drawnow

end






function [cone_weights,cone_colors] = find_best_weights(rfs,cone_info,kernel_params)

% more convenient variable
cone_centers = cone_info(:,end-2:end-1);

num_cones = size(cone_info,1);
cone_fit_size = [size(rfs,1) size(rfs,2)];


% get bw kernel for each cone
bw_kernel = zeros(cone_fit_size(1),cone_fit_size(2),num_cones);
for cc = 1:num_cones
    % make bw kernel
    bw_kernel(:,:,cc) = make_gaussian(struct('center',cone_centers(cc,1:2),...
        'y_size',cone_fit_size(1),'x_size',cone_fit_size(2),...
        'center_radius',kernel_params.center_radius,'normalize','sum'));
end



%get best color for each cone

% identify each cone's RGB
for cc = 1:num_cones
    summed_rgb(cc,:) = sum(reshape(reshape(rfs,[],size(rfs,3)*size(rfs,4))'*...
        reshape(bw_kernel(:,:,cc),[],1),size(rfs,3),size(rfs,4)),2);
end

% identify which cone types those are
cone_types = classify_cones(summed_rgb, kernel_params.cone_rgb, 'algorithm','nearest line');

% 0 for L/M, 1 for S
cone_colors = (cone_types=='S');



% compute matrix to translate cone weights to RF fit
% # weights = (# cones) * (# RFs)
% # dims in RFs = (# RFs) * (# colors) * X * Y
M = [];

% add in values for each cone
for cc = 1:num_cones

    % add color
    rgb_kernel = zeros(size(bw_kernel,1),size(bw_kernel,2),size(rfs,3));
    % get rgb weights for this cone
    rgb = cone_rgb(cone_colors(cc),kernel_params.cone_rgb);
    %rgb = summed_rgb(cc,:);
    % fill in color
    for dd = 1:size(rfs,3)
        % weight cone kernel for this color
        rgb_kernel(:,:,dd) = bw_kernel(:,:,cc)*rgb(dd);
    end

    % set norm to 1
    rgb_kernel = rgb_kernel/(norm(reshape(rgb_kernel,[],1)).^0.5);

    % add to accumulating matrix
    clear temp
    for rr=1:size(rfs,4)
        temp{rr}=reshape(rgb_kernel,[],1);
    end

    M = [M blkdiag(temp{:})];
end


% compute weights for each cone
%cone_weights = M'*reshape(rfs,[],1);
cone_weights = M\reshape(rfs,[],1);

%rf_fits = M*cone_weights;
rf_fits = reshape(M*cone_weights,size(rfs));

cone_weights = reshape(cone_weights,size(rfs,4),num_cones)';




% plot best fitting weights
if 1

    figure(14);clf
    plot_axes = subplot_axes(14,[0 0 1 1],0.1,0.1,size(rfs,4),1);

    disp(sum((reshape(rfs-rf_fits,[],1)).^2))

    % plot them
    pa = plot_axes;
    for rr = 1:size(rfs,4)

        % plot RFs
        cla(pa{rr})
        %imagesc(norm_image([rfs(:,:,:,rr);rf_fits(:,:,:,rr);rf_fits_p1(:,:,:,rr);rf_fits_p2(:,:,:,rr)]),'Parent',pa{rr});
        imagesc(norm_image([rfs(:,:,:,rr);rf_fits(:,:,:,rr)]),'Parent',pa{rr});axis(pa{rr},'image')
        % plot cones
        hold(pa{rr},'on')
        % over RFs
        plot(pa{rr},cone_info(:,end-2),cone_info(:,end-1),'.k')
        % over fits
        plot(pa{rr},cone_info(:,end-2),cone_info(:,end-1)+size(rfs,1),'.k')
    end
    drawnow

end

