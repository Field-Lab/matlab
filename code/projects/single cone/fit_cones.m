function [fit_cone_info,rf_fit] = fit_cones(rf,cone_info,varargin)
% FIT_CONES     fit ideal cone shapes to cones in a RF
%
% usage:  [fit_cone_info,rf_fit] = fit_cones(rf,cone_locations,params)
%
% arguments:        rf - 2D intensity map of RF
%            cone_info - Nx3 matrix, initial guesses of cone locations
%                           first column x coord, second y, third strength
%
% output:   fit_cone_info - Nx3 matrix, fits to cone locations, some format as cone_info
%                  rf_fit - final fit to the RF
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
%                                           center_radius       1       center gaussian radius
%                                           surround_radius     2.5     surround gaussian radius
%                                           surround_scale      0       surround gaussian scale (NOT relative)
%                                           effective_radius    2       region in which to compute the gaussian
%   
%                                       
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
p.addParamValue('algorithm', 'lsqcurvefit',@(x)any(strcmpi(x,{'fminsearch','lsqcurvefit','com'})));

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
            defaults.center_radius = 1.5;
            defaults.surround_radius = 2.5;
            defaults.surround_scale = 0;
            defaults.effective_radius = 2;
            
        otherwise
            error('Cone kernel type ''%s'' not valid',params.cone_kernel.type)
    end

    % combine user and default parameters
    params.cone_kernel = default_params( defaults, params.cone_kernel);
end




% set up figure
if ~isempty(params.figure)
    if params.figure > 0
        figure(params.figure)
    else
        params.figure = figure;
    end
    figure(params.figure);clf
end


% % optimize fit, save cone info for best fit
% fit_cone_info = fminsearch(@(cone_info)get_fit_err(rf,cone_info,params.cone_kernel,params.figure),...
%     cone_info,... initial parameters
%     optimset(params.optim{:})); % parameters specifying tolerance threshold

% optimize fit, save cone info for best fit
switch params.algorithm
    case 'fminsearch'
        fit_cone_info = fminsearch(@(cone_info)get_fit_err(rf,cone_info,params.cone_kernel,params.figure),...
            cone_info,... initial parameters
            optimset(params.optim{:})); % parameters specifying tolerance threshold
    case 'lsqcurvefit'
        % make function that reads in the parameters AND the x data
        fit_fcn = @(cone_info,xx)get_rf_fit(rf,cone_info,params.cone_kernel,params.figure);

%         size(1)
%         size(rf)
%         pause
        
        fit_cone_info = lsqcurvefit(...
            fit_fcn,...
            cone_info,... initial parameters
            1,... x data (not relevant here)
            rf,... y data
            [],[],optimset(params.optim{:})); % y data
    otherwise
        error('fitting algorithm ''%s'' not recognized.',params.algorithm)
end

% return the fit with the best parameters
rf_fit = get_rf_fit(rf,fit_cone_info,params.cone_kernel,params.fig_final);





% take in cone locations and return error rate 
function err = get_fit_err(rf,cone_info,kernel_params,fig)

% get RF fit based on current cone info
rf_fit = get_rf_fit(rf,cone_info,kernel_params,fig);

% compute rmse
err = sqrt(mean(mean((rf_fit - rf).^2)));




% return a fit
function rf_fit = get_rf_fit(rf,cone_info,kernel_params,fig)

% make matrix to store fit
rf_fit = zeros(size(rf));

cone_fit_size = size(rf);

% compute sum of cone fits
for cc = 1:size(cone_info,1)
    
    % generate cone fit
    cone_fit = make_gaussian(struct('center',cone_info(cc,1:2),...
        'y_size',cone_fit_size(1),'x_size',cone_fit_size(2),...
        'center_radius',kernel_params.center_radius,'normalize','sum'));
    
    % add to rf fit
    rf_fit = rf_fit + cone_info(cc,3)*cone_fit;
    
end

% plot current fit
if ~isempty(fig)
    switch 3
        case 1
            figure(fig);clf;colormap('gray');
            imagesc(rf);hold on; plot(cone_info(:,1),cone_info(:,2),'.r');drawnow
        case 2
            figure(fig);clf;colormap('gray');
            subplot(121);imagesc(rf); hold on; plot(cone_info(:,1),cone_info(:,2),'.r')
            subplot(122);imagesc(rf_fit); hold on; plot(cone_info(:,1),cone_info(:,2),'.r')
            %drawnow
        case 3
            figure(fig);clf;colormap('gray');
            imagesc([rf rf_fit]);hold on; plot(cone_info(:,1),cone_info(:,2),'.r');drawnow
            
    end
end

