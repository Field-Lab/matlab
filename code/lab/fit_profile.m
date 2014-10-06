function [fit_params, the_fit] = fit_profile(x,y, varargin)
% fit_profile     fit an analytic shape to a set of x and y values
%
% usage:  [fit_params, the_fit] = fit_profile(x,y, varargin)
%
% arguments:      x,y - x and y values to fit
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    fit_params - struct of parameters
%                the_fit - shape of the fit, same length as x (should be similar to y)
%
%
% optional params, their default values, and what they specify:
%
% verbose           false               show output
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure. if empty, don't plot
% type              'dog'               which shape to fit
%                                           'dog' - difference of gaussians
% algorithm         'fminsearch'        fitting algorithm
%                                           'fminsearch'
%                                           'lsqcurvefit'
% optim             {'TolX',0.1,'Display','off'}
%                                       arguments to optimset to control optimization of fit.
%                                           see documentation of fminsearch for options.
%
% additional params for type == 'dog'
%
%   center_radius           see make_gaussian for explanation of parameters
%   center_scale            
%   surround_radius
%   surround_scale
%   center
%   fit_center_radius       allow center_radius to vary in the fit
%   fit_center_scale        etc...
%   fit_surround_radius
%   fit_surround_scale
%   fit_center
%
% NOTE: dog fits are subject to these constraints:
%
%           center_radius >= 0
%           surround_scale >= 0
%           surround_radius >= center_radius
%
%
%
% gauthier 2008-10
%




% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional arguments
p.addParamValue('verbose', false);
p.addParamValue('fig_or_axes', []);
p.addParamValue('type','dog', @(x)any(strcmpi(x,{'dog'})));
p.addParamValue('optim',{'TolX',0.1,'Display','off'});
%p.addParamValue('algorithm', 'fminsearch',@(x)any(strcmpi(x,{'fminsearch','lsqcurvefit'})));
p.addParamValue('algorithm', 'lsqcurvefit',@(x)any(strcmpi(x,{'fminsearch','lsqcurvefit'})));

% forked parameters
%    type = dog
p.addParamValue('center_radius','default value');
p.addParamValue('center_scale','default value');
p.addParamValue('surround_radius','default value');
p.addParamValue('surround_scale','default value');
p.addParamValue('center','default value');
p.addParamValue('fit_center_radius','default value');
p.addParamValue('fit_center_scale','default value');
p.addParamValue('fit_surround_radius','default value');
p.addParamValue('fit_surround_scale','default value');
p.addParamValue('fit_center','default value');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% set parameters forked from type

% set default values in an inputParser object
p_temp = inputParser;
switch params.type
    case 'dog'
        % by default, fit one gaussian with fixed center
        p_temp.addParamValue('center_radius',1);
        p_temp.addParamValue('center_scale',[]); % if empty, will make a guess based on the maximum value
        p_temp.addParamValue('surround_radius',1);
        p_temp.addParamValue('surround_scale',0);
        p_temp.addParamValue('center',0);
        p_temp.addParamValue('fit_center_radius',true,@islogical);
        p_temp.addParamValue('fit_center_scale',true,@islogical);
        p_temp.addParamValue('fit_surround_radius',false,@islogical);
        p_temp.addParamValue('fit_surround_scale',false,@islogical);
        p_temp.addParamValue('fit_center',false,@islogical);
end

% add forked parameters to the params struct
params = add_forked_params(params,p_temp);


% BODY OF THE FUNCTION

% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% show output
if params.verbose
    start_time = clock; % note when it started
end


    

% ensure x and y are row vectors
x = reshape(x,1,[]);
y = reshape(y,1,[]);


% make initial guess about center scale
if isempty(params.center_scale)
    params.center_scale = max(y);
end


% get initial conditions for fit
input_params = nums_from_struct(params);

% note which to fit
which_to_fit = [...
    params.fit_center_radius...
    params.fit_center_scale...
    params.fit_surround_radius...
    params.fit_surround_scale...
    params.fit_center ];

% note which indices in the vector are fit, and which are not
fit_indices = find(which_to_fit);
fixed_indices = find(~which_to_fit);


% optimize parameters which are to be fit, locking in initial values for parameters which are not fit
switch params.algorithm
    case 'fminsearch'
        best_fit_params = fminsearch(...
            @(fit_params)get_fit_err(x,y,...
            fit_params,input_params(fixed_indices),fit_indices,fixed_indices,plot_axes),...
            input_params(fit_indices),... initial parameters
            optimset(params.optim{:})); % parameters specifying tolerance threshold
    case 'lsqcurvefit'
        % make function that reads in the parameters AND the x data
        fit_fcn = @(fit_params,xx)get_fit_only(xx,y,...
            fit_params,input_params(fixed_indices),fit_indices,fixed_indices,plot_axes);
        
        best_fit_params = lsqcurvefit(...
            fit_fcn,...
            input_params(fit_indices),... initial parameters
            x,...x  data
            y,...
            [],[],optimset(params.optim{:})); % y data
    otherwise
        error('fitting algorithm ''%s'' not recognized.',params.algorithm)
end


% combine locked parameters and fit parameters
final_params = input_params;
final_params(fit_indices) = best_fit_params;


% compute the fit with the best parameters
fit_x = x;
fit_y = make_gaussian('dim',1,'x_vals',fit_x,struct_from_nums(final_params));

% make variables to return
the_fit = [fit_x; fit_y]';
fit_params = struct_from_nums(apply_constraints(final_params));



% display how long it took
if params.verbose
    fprintf('fit %d parameters on %d points in %0.1f seconds\n',...
        length(fit_indices),length(x),etime(clock,start_time));
end





% function that takes in fit parameters and returns the fit
function the_fit = get_fit_only(x,y,fit_params,fixed_params,fit_indices,fixed_indices,plot_axes)

% combine fit and fixed params
all_params(fit_indices) = fit_params;
all_params(fixed_indices) = fixed_params;

% apply constraints
all_params = apply_constraints(all_params);

% get RF fit based on current parameters
the_fit = make_gaussian('dim',1,'x_vals',x,struct_from_nums(all_params));

% plot the fit
if ~isempty(plot_axes)
    % sort points
    [x,ii] = sort(x);
    plot(plot_axes,x,y(ii),'.r',x,the_fit(ii),'k');
    drawnow
end





% function that takes in fit parameters and returns error rate
function err = get_fit_err(x,y,fit_params,fixed_params,fit_indices,fixed_indices,plot_axes)

% get the fit
the_fit = get_fit_only(x,y,fit_params,fixed_params,fit_indices,fixed_indices,plot_axes);

% compute rmse
err = sqrt(mean(mean((y - the_fit).^2)));




function fp = struct_from_nums(fit_params)
% put list of fit params into a struct

fp.center_radius = fit_params(1);
fp.center_scale = fit_params(2);
fp.surround_radius = fit_params(3);
fp.surround_scale = fit_params(4);
fp.center = fit_params(5);



function nums = nums_from_struct(params)
% put struct of fit params into a list

nums = [...
    params.center_radius ...
    params.center_scale...
    params.surround_radius...
    params.surround_scale...
    params.center];



function all_params = apply_constraints(all_params)

% ensure center_radius >= 0
if all_params(1) < 0; all_params(1) = realmin; end

% ensure surround_radius >= center_radius
if all_params(3) < all_params(1); all_params(3) = all_params(1); end

% ensure surround_scale >= 0
if all_params(4) < 0; all_params(4) = 0; end




