function [sig_stixels,rf_strength_out] = find_cones_in_rf(rf, varargin)
% FIND_CONES     find cones in a RF
%
%   1. make noise sigma = 1
%   2. apply spatial filter (optional)
%   3. compute RF strength at each stixel (collapse across RGB)
%   4. return stixels which seem not to come from the noise distribution
%   
%
% usage:  [sig_stixels,rf_strength_out] = find_cones_in_rf(rf)
%
% arguments:  rf - 2D or 3D RF (x,y,color)
%
% output:   sig_stixels - logical matrix same x-y size as rf, showing which stixels are significantly sensitive
%           rf_strength_out - scalar matrix same x-y size as rf, showing RF sensitivity at each point
%
% optional fields in params, their default values, and what they specify:
%
% fig_rfs           []          which figure to plot in.  if 0, make new.  if empty, don't plot
%
% strength          {'inner',[.3 1 .3]}
%                               how to get the RF strength at each stixel
%                                   'vector length' - square each color channel, sum them, and take square root
%                                   {'inner',v}     - compute inner product with length-3 vector v
%                               
%
% filter            struct('type','given','filt',make_gaussian)
%                               struct describing the spatial filter to apply before finding significant stixels
%                                   if empty, no filter is applied
%                                   the spatial filter is described by the fields of a struct
%                                   only required field is 'type'
%                                   all other fields are optional, with different default values depending on the type of filter
%
%                                   for type = 'dog', these are the optional fields and their default values
%                                           center_radius       1       center gaussian radius
%                                           surround_radius     2.5     surround gaussian radius
%                                           surround_scale      0.3     surround gaussian scale relative to center
%                                           x_size              5       size of filter matrix in x
%                                           y_size              5       in y
%                                           center              [3 3]   location of center
%
%
% selection         struct('type','max')
%
%                               struct describing which algorithm to use to select significant stixels
%                                   only required field is 'type', and it may take on these values:
%                                   
%                                   'thresh' - apply simple threshold
%                                   'max'    - only stixels that are at a local maximum in their neighborhood
%                                   'lasso'  - reconstruct RF from unitary cones, while limiting L1 norm
%                                   
%                                   all other fields are optional, with different default values depending on the algorithm
%                                   
%                                   for type = 'thresh', these are the optional fields and their default values
%                                           thresh  	5           threshold, see below
%                                   
%                                   for type = 'max'
%                                           thresh  	5           threshold, see below
%                                           radius  	1           radius of search window in which to look for local maxima
%                                           maxima  	'center'    how to select maxima
%                                                                       'center' - must be in center of window
%                                                                       'any'    - can be anywhere in search window
%                                   
%                                   for type = 'lasso'
%                                                     
%
%
%
%
%   note on computing threshold
%
%   threshold is specified in terms of SNR.  The notion of SNR assumes a gaussian noise distribution.  This is
%   generalized to other distributions by computing the probability that selected pixels are not from the noise distribution.
%
%   for a Gaussian distribution,
%       sigma is estimated and the absolute threshold is T * sigma  (T = user specified threshold value)
%   for a chi-square distribution,
%       threshold is computed in 2 steps.  First, the p value is computed for the specified threshold.
%       Second, this p value is translated to the chi square value.  Thus the absolute threshold is chi2inv(normcdf(T),df)
%   for an unknown distribution,
%       it is assumed to be approximately Gaussian, and the "left half" (values < mean) is used to estimate sigma
%       In practice, this works pretty well.
%
%
% gauthier 2008-09
%



%   OLD NOTES.  IGNORE FOR NOW (UNTIL WE CHANGE IT BACK TO THIS WAY!)
%
%   note on computing threshold
%
%   threshold is defined in terms of the probability that selected pixels are not from the noise distribution.
%   for a Gaussian distribution,
%       sigma is estimated and the threshold is computed in terms of sigmas (norminv)
%   for a chi-square distribution,
%       threshold is computed using the inverse cumulative chi square function (chi2inv)
%   for an unknown distribution,
%       it is assumed to be approximately Gaussian, and the "left half" (values < mean) is used to estimate sigma
%       needless to say, this is not ideal, but it works pretty well
%
%   to use perhaps more familiar units of "snr", use normcdf(SNR) as the argument
%
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('fig_rfs', []);
p.addParamValue('strength', 'vector length');
%p.addParamValue('strength', {'inner',[.3 1 .3]});
%p.addParamValue('strength', {'inner or',[0.4044 0.8854 0.2292;0.1483 0.9161 0.3726;0.0174 0.0792 0.9967]});
p.addParamValue('filter',struct('type','given','filt',make_gaussian));
p.addParamValue('selection',struct('type','max'));
p.addParamValue('plot_filter', false);

% parameters to be passed on
%    space summary function
% p.addParamValue('space_bins', 'default value');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% generate structs to pass on
% time_params = make_struct_to_pass(p.Results,{'color','time_color','bins',{'time_bins','bins'}});

        
% get default filter parameters
clear defaults
if ~isempty(params.filter)
    switch params.filter.type

        case 'dog'
            defaults.type = 'dog';
            defaults.center_radius = 1.0;
            defaults.surround_radius = 2.5;
            defaults.surround_scale = 0.3;
            defaults.x_size = 5;
            defaults.y_size = 5;
            defaults.center = [1 1]*3;
            
        case 'given'
            defaults.type = 'given';
            defaults.filt = make_gaussian;

        otherwise
            error('Filter type ''%s'' not valid',params.filter.type)
    end
    
    % combine user and default parameters ("fp" for "filter parameters")
    fp = default_params( defaults, params.filter);
end


% get default selection parameters
clear defaults
switch params.selection.type

    case 'thresh'
        defaults.type = 'thresh';
        defaults.thresh = 5;

    case 'max'
        defaults.type = 'max';
        defaults.thresh = 5;
        defaults.radius = 1;
        defaults.maxima = 'center';

    case 'lasso'
        defaults.type = 'lasso';

    otherwise
        error('Selection type ''%s'' not valid',params.selection.type)
end

% combine user and default parameters (sp = selection parameters)
sp = default_params( defaults, params.selection);





% % set parameters forked from dim 
% 
% % set default values in an inputParser object
% p_temp = inputParser;
% switch p.Results.dim
%     
%     case 2 % 2D gaussian
%         p_temp.addParamValue('x_size',3);
%         p_temp.addParamValue('y_size',3);
%         p_temp.addParamValue('center',[2 2],@(x)numel(x) == 2);
%         p_temp.addParamValue('effective_radius', Inf);
%         
%     case 1 % 1D gaussian
%         p_temp.addParamValue('x_vals',-5:0.1:5);
%         p_temp.addParamValue('center',0,@(x)numel(x) == 1);
% end
% 
% % add forked parameters to the params struct
% params = add_forked_params(params,p_temp);




% set up figure
if ~isempty(params.fig_rfs)
    if params.fig_rfs > 0
        figure(params.fig_rfs)
    else
        params.fig_rfs = figure;
    end
    sub_x = 3;
    sub_y = 2;
    figure(params.fig_rfs);clf
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZE NOISE SIGMA TO 1 

rf_noise_sigma = robust_std(reshape(rf,[],1));
if rf_noise_sigma > 0
    rf = rf/rf_noise_sigma;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLY SPATIAL FILTER

if ~isempty(params.filter)

    % get specified filter
    switch params.filter.type
        case 'dog'
            
            % generate DOG filter

            % compute distance of each point from the center
            x_dist = ([1:fp.x_size] - fp.center(1)).^2;
            y_dist = ([1:fp.y_size] - fp.center(2)).^2;
            dist = sqrt(ones(fp.y_size,1)*x_dist + y_dist'*ones(1,fp.x_size));

            % generate profile function for center and surround
            filt_fcn = @(x)( exp(-1*(x/fp.center_radius).^2) - fp.surround_scale*exp(-1*(x/fp.surround_radius).^2) );

            % apply profile to these distances
            filt = filt_fcn(dist);
            
        case 'given'
            
            filt = params.filter.filt;
            
                
        otherwise
            error('Filter type not recognized: %s',params.filter.type)
    end
    
    % normalize filter so that sum of squares = 1
    filt = filt ./ sqrt(sum(sum(filt.^2)));

    % apply filter to RF vector length in each color
    rf_filt = zeros(size(rf));
    for cc=1:size(rf,3)
        rf_filt(:,:,cc) = filter2(filt,rf(:,:,cc));
    end

    
    % plot filter
    if ~isempty(params.fig_rfs) && params.plot_filter
        figure(params.fig_rfs+1);clf;imagesc(norm_image( filt ));
    end
    

else % if no filter

    % make the variable "rf_filt" equal to the rf vector length
    rf_filt = rf;

end


% get number of colors
num_col = size(rf_filt,3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE RF STRENGTH
%
% also define function that converts threshold to relevant units
%
% rf_strength contains values that are output
% rf_test contains values which will be tested
% f_thresh is the function applied to the threshold to determine which points to keep


% determine how to compute strength
if ischar(params.strength)
    strength_type = params.strength;
elseif iscell(params.strength)
    strength_type = params.strength{1};
end

switch strength_type
    case 'vector length'

        % get vector length of each pixel
        rf_strength = sqrt(sum(rf_filt.^2,3));
        % then the SQUARE of noise pixels will follow a chi square distribution (d.f. = # color dimensions of the original rf)
        
        % rf test = rf strength
        rf_test = rf_strength;

        % base threshold on chi square distribution
        f_thresh = @(th)sqrt(chi2inv(normcdf(th),size(rf,3)));

        
    case 'inner'
        
        % put into convenient variable
        v = params.strength{2};
        
        % ensure length of v = number of colors
        if length(v) ~= num_col
            error('Vector for dot product has length (%d) not equal to the number of colors (%d).',length(v),num_col)
        end
        
        % apply dot product to each stixel
        
        % reshape so rf has one color per column
        rf_temp = reshape(rf_filt,[],num_col);
        
        % take dot product
        rf_dot_prod = rf_temp*v';
        
        % normalize to make sigma 1
        % if the different color channels are independent, this is unnecessary
        % but if they are correlated, then sigma will not be 1
        rf_dot_prod = rf_dot_prod/robust_std(rf_dot_prod);
        
        % reshape
        rf_strength = reshape(rf_dot_prod,size(rf,1),size(rf,2));
        
        % normalize so sigma is 1
        rf_strength = rf_strength / sqrt(sum(v.^2));
        
        % take absolute value
        rf_strength = abs(rf_strength);
        
        % rf test = rf strength
        rf_test = rf_strength;
        
        % account for the absolute value
        f_thresh = @(th)(-norminv(normcdf(-th)/2));
        
    case 'inner or'
        % rf_test will be boolean
        % f_thresh will return 0.5
        % rf_strength will be the maximum across the various inner products
        
        % get vectors for dot product
        vectors = params.strength{2};
        
        % normalize
        vectors = normalize_to_unit_sphere(vectors);
        
        % reshape so rf has one color per column
        rf_temp = reshape(rf_filt,[],num_col);
        
        % go through each vector
        for vv =1:size(vectors,1)
            
            % compute dot product
            temp = rf_temp*vectors(vv,:)';
            
            % normalize to make sigma 1
            % if the different color channels are independent, this is unnecessary
            % but if they are correlated, then sigma will not be 1
            temp = temp/robust_std(temp);
            
            % take absolute value
            rf_dot_prod_temp(:,:,vv) = abs(reshape(temp,size(rf,1),size(rf,2)));
            
            % note whether absolute value exceeds threshold
            rf_tests_temp(:,:,vv) = rf_dot_prod_temp(:,:,vv) > -norminv(normcdf(-sp.thresh)/2.2/2);
            
            % threshold is adjusted in a complicated way
            % / 2 to account for absolute value
            % / 2.2 because, empircally, this is required to give the appropriate false positive rate
            
        end
        
        rf_test = logical(sum(rf_tests_temp,3));
        
        rf_strength = max(rf_dot_prod_temp,[],3);
        
        f_thresh = @(th)0.5;
        
    otherwise
        error('Method to get RF strength ''%s'' not recognized.',strength_type)
        
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET SIGNIFICANT STIXELS


% define threshold for params that use it
switch params.selection.type 
    case {'thresh','max'}
        abs_thresh = f_thresh(sp.thresh);
end


% look for significant stixels
cone_kernel = []; % default
switch params.selection.type

    case 'thresh' % apply simple threshold

        % return significant stixels
        sig_stixels = rf_test > abs_thresh;
        
        % return rf strength as vector length
        rf_strength_out = rf_strength;


    case 'max'  % return any pixel that is a maximum in its local neighborhood
        
        
        sig_stixels = find_local_maxima(rf_test,'thresh',abs_thresh,'radius',sp.radius,'return','matrix');
        
        % return rf strength as vector length
        rf_strength_out = rf_strength;
        
        
        
    case 'lasso'
        
        % choose region of interest
        
        switch 1
            case 1  % simple radius around COM
                roi_rad = 20;
                ctr = round(rf_com(rf,struct('ss_params',struct('select','thresh','thresh',4.5))));
                roi_x = ctr(1)+[-roi_rad:roi_rad];
                roi_y = ctr(2)+[-roi_rad:roi_rad];
                
            case 2  % filter with gaussian, look for signal
                
        end
        
        % rf_strength follows a sqrt(chi_sqr) distribution.  rf_filt should be similar.
        rf_filt = rf_filt .^2;
        % re-zero rf_filt
        rf_filt = rf_filt - median(reshape(rf_filt,[],1));
        
        % take RF in this region
        %rf_roi = rf_filt(roi_y,roi_x);
        if length(size(rf)) == 3
            rf_roi = rf(roi_y,roi_x,2);
        elseif length(size(rf)) == 2
            rf_roi = rf(roi_y,roi_x);
        else
            error('size of rf not supported')
        end
        if 1
            figure(params.fig_rfs+5)
            hist(reshape(rf_filt(:,1:100),[],1),50)
            drawnow
        end
        
        % compute cone kernels
        %[W_cone_kernels,cone_kernel] = compute_cone_kernels(size(rf_roi,2),size(rf_roi,1));
        [W_cone_kernels,cone_kernel] = compute_cone_kernels(size(rf_roi,2),size(rf_roi,1),...
            struct('figure',5,'kernel_params',struct('type','dog','center_radius',0.8)));
        
        % get expected strength of RF by adding all stixels
        expected_str = sum(sum(rf_roi));
        
        % divide by cone kernel weights to get t parameter
        t = expected_str / sum(sum(cone_kernel));
        
        % and make it smaller, just to be on the safe side
        t = 1.5*t;
        
        % compute lasso
        lasso_reconstruction = LassoActiveSet(W_cone_kernels,reshape(rf_roi,[],1),t);
        
        % embed within larger matrix
        rf_strength_out = zeros(size(rf_filt));
        rf_strength_out(roi_y,roi_x) = reshape(lasso_reconstruction,length(roi_y),length(roi_x));
        
        % set signficant stixels
        sig_stixels = (rf_strength_out > 0);
            
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THRESHOLD



                    
% % compute threshold based on noise sigma, assuming "left half" is all noise
% 
% % get "left half" distribution
% rf_filt_vals = reshape(rf_filt,[],1);
% less_than_mean = rf_filt_vals(rf_filt_vals < mean(rf_filt_vals)) - mean(rf_filt_vals);
% % get sigma of these points (include a matching "right half")
% noise_sigma = robust_std([less_than_mean' -less_than_mean']);
% 
% % compute how many sigmas is the equivalent of the specified threshold
% num_sigmas = sp.thresh;
% 
% % take threshold as mean plus a number of noise sigmas
% abs_thresh = mean(rf_filt_vals) + num_sigmas * noise_sigma;





% optionally plot rf before and after filtering

if params.fig_rfs

    % set up figure for plotting cell
    figure(params.fig_rfs)

    % plot original RF
    a1 = subplot(sub_y,sub_x,1);
    imagesc(norm_image(rf)); axis image
    title('original')

    % vector length
    a2 = subplot(sub_y,sub_x,2);
    imagesc(norm_image(rf_strength)); axis image
    title('vector length')

    % filtered version
    a3 = subplot(sub_y,sub_x,3);
    imagesc(norm_image(rf_filt)); axis image
    title('filter')

    % filtered version, points exceeding threshold
    a4 = subplot(sub_y,sub_x,6);
    imagesc(norm_image(rf_filt > abs_thresh)); axis image
    title('filter, thresh')

    % UNfiltered version, points exceeding threshold
    a5 = subplot(sub_y,sub_x,5);
    imagesc(norm_image(rf_strength > abs_thresh)); axis image
    title('thresh')

    axes_to_link_1 = [a1 a2 a3 a4 a5];

end



% optionally plot significant stixels

if params.fig_rfs

    % go to figure for plotting cell
    figure(params.fig_rfs)

    % plot significant stixels
    a4 = subplot(sub_y,sub_x,4);
    imagesc(norm_image(sig_stixels)); axis image
    title('significant stixels')
    axes_to_link_1 = [axes_to_link_1 a4];

    
    % plot lasso stuff
    if strcmp(params.selection.type,'lasso')

        % plot lasso cone weights
        a5 = subplot(sub_y,sub_x,5);
        imagesc(norm_image(rf_strength_out)); axis image
        title('lasso cone weights')
        axes_to_link_1 = [axes_to_link_1 a5];

        % plot lasso reconstrutcion
        a6 = subplot(sub_y,sub_x,6);
        imagesc(norm_image(imfilter(rf_strength_out,cone_kernel))); axis image
        title('lasso reconstruction')
        axes_to_link_1 = [axes_to_link_1 a6];
    end

        
    
    % zoom in on the cell
    ss_params.thresh = 3.5;
    rf_ctr = rf_com(rf,'ss_params', ss_params);
    zoom_rad = 35;
    set(gca,'XLim',[rf_ctr(1)-zoom_rad rf_ctr(1)+zoom_rad],'YLim',[rf_ctr(2)-zoom_rad rf_ctr(2)+zoom_rad ])
    linkaxes(fliplr(axes_to_link_1));

end


% optionally plot values and threshold
if 0 && ~isempty(params.fig_rfs)

    % make figure
    figure(params.fig_rfs+3);clf

    % plot square of vector length
    hist(reshape(rf_strength(:,1:100).^2,[],1),200)

    % plot square of threshold
    hold on; plot([1 1]*abs_thresh.^2,get(gca,'YLim'),'k')

    % plot ideal chi square
    xx = 0:.2:50;
    plot(xx,chi2pdf(xx,size(rf,3))*max(get(gca,'YLim'))*4,'r')

    % make figure showing if original rf was Gaussian
    %figure(params.fig_rfs+4);clf

    %normplot(reshape(rf(:,1:180,3),[],1))

end


