function [sig_stixels, params, rf_strength_out] = significant_stixels(sta, varargin)
% significant_stixels     identify "significant" stixels in an STA or RF 
%
% strategy:
%
%  * spatially filter the STA (if desired)
%  * collapse across time
%  * normalize so that in each color channel mean = 0 and sigma = 1
%  * reduce to 1 color
%  * select stixels that exceed threshold
%
%
%
% usage:  [sig_stixels, params, rf_strength] = significant_stixels(sta, params)
%
% arguments:      sta - 2d, 3d, or 4d matrix (height, width, color, time)
%              params - struct of optional parameters (see below)
%
% outputs:             sig_stixels - YxX logical matrix, stixels which are "significant"
%                           params - struct containing the parameters that were used
%                      rf_strength - the value at each stixel which determined whether it was significant
%                                       note that these values must be compared to the absolute threshold to make sense
%
%
% optional fields in params, their default values, and what they specify:
%
% filter      	[]          filter to convolve with each color channel in each frame
% time          'std'       how to collapse across time in each color channel
%                               'std' - take std
%                               'max' - take maximum deviation from 0
%
% strength      {'inner',[1 1 1]}
%                           how to combine information across color channels.
%                           can be a string, or a cell array with the first item a string,
%                           with remaining parameters following:
%                               'vector length' - sqrt of sum of squares
%                               {'inner',v}     - take inner product with vector v
%                               {'inner or',V}  - take inner product with vectors V(i,:)
%                                                   keep a stixel if it exceeds threshold for ANY i
% select        'thresh'    how to choose which stixels are significant
%                               'thresh' - apply simple threshold
%                               'max'    - only keep local maxima
%
% robust_std_method 3       which algorithm to use for calculating
%                               robust_std; see robust_std for more
%                               information.
% 
% these arguments only apply if select == 'thresh'
%       thresh      5       threshold in "SNR units".  actual threshold is converted to appropriate units
%                               based on the equivalent p value.
% 
% these arguments only apply if select == 'max'
%       thresh      5       same units as above.  This threshold is applied prior to looking for local maxima.
%       radius      1       radius in which the stixel must be the maximum
%                               see find_local_maxima for details.
%                       
%
%
% 2008-10 gauthier    
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('filter', []);
p.addParamValue('time', 'std',@(x)any(strcmpi(x,{'max','std'})));
p.addParamValue('select','thresh',@(x)any(strcmpi(x,{'thresh','max'})));
p.addParamValue('robust_std_method', 3);

% choose default 'strength' value based on number of colors
if size(sta,3) == 1
    p.addParamValue('strength', {'inner',1});
else
    %p.addParamValue('strength', {'inner',[.3 1 .3]});
    p.addParamValue('strength', {'inner',[1 1 1]});
    %p.addParamValue('strength', {'inner or',[0.4044 0.8854 0.2292;0.1483 0.9161 0.3726;0.0174 0.0792 0.9967]});
end

%    selection  == 'thresh' or 'max'
p.addParamValue('thresh','default value');
%    selection  == 'max'
p.addParamValue('radius','default value');



% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% set parameters forked from selection type 

% set default values in an inputParser object
p_temp = inputParser;
switch params.select
    
    case 'thresh'
        p_temp.addParamValue('thresh',5);
        
    case 'max'
        p_temp.addParamValue('thresh',5);
        p_temp.addParamValue('radius',1);
end

% add forked parameters to the params struct
params = add_forked_params(params,p_temp);




% if STA is empty, return empties
if isempty(sta)
    sig_stixels = [];
    params = struct;
    rf_strength_out = [];
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLY SPATIAL FILTER

if ~isempty(params.filter)

    % put into more convenient variable
    filt = params.filter;
    
    % normalize filter so that sum of squares = 1
    filt = filt ./ sqrt(sum(sum(filt.^2)));

    % apply filter to RF vector length in each color
    sta_filt = imfilter(sta,filt);

else % if no filter

    % make the variable "sta_filt" equal to the sta
    sta_filt = sta;

end


% get number of colors
num_col = size(sta_filt,3);










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COLLAPSE ACROSS TIME
%
% only if needed

% figure(41)
% imagesc(norm_image(sta(:,:,:,28)))
% drawnow

if size(sta,4) > 1
    
    % normalize values
    %
    % get mean 0, std 1 in each color channel

    % reshape to have each stixel's colors in a row
    sta_reshape = reshape(permute(sta_filt,[1 2 4 3]),[],size(sta,3));

    % recenter each color's distribution at 0
    %sta_reshape = sta_reshape - repmat( median(sta_reshape) ,size(sta_reshape,1),1);
    
    % get noise sigma of each color channel
    %figure
    %hist(sta_reshape(:,2))
    noise_sigmas = robust_std(sta_reshape, params.robust_std_method);
    
    % introduced by GDF to handle case where robust_std returns zeros
%     if any(noise_sigmas == 0)
%         % find peak value of STA
%         max_val = max(abs(sta_reshap));
%         % make noise (sigma) level 0.05% of this peak value
%         noise_sigmas(1:size(sta_3)) = max_val ./ 0.05
%         pause
%     end
        

    % divide by noise sigmas
    sta_reshape = sta_reshape./repmat(noise_sigmas,size(sta_reshape,1),1);

    % reshape
    sta_renorm = permute(reshape(sta_reshape,size(sta,1),size(sta,2),size(sta,4),size(sta,3)),[1 2 4 3]);



    switch params.time
        case 'max'
            % get the largest deviation from zero at each color/stixel
            rf_highs = max(sta_renorm,[],4);
            rf_lows = max(-sta_renorm,[],4);
            rf = max(rf_highs,rf_lows).*sign(rf_highs-rf_lows);
            
        case 'std'
            % take the standard deviation across time
            rf = std(sta_renorm,1,4);
            
            % transform to normal distribution
            
            % the distribution of standard deviation values with 1/n normalization (s_hat) has the following property:
            %   n*(s_hat^2)/(s_true^2)  =  chi square with n-1 degrees of freedom
            %       where n is the number of samples and s_true is the true standard deviation
            % here s_true is 1 (by the normalization above), and n = size(sta,4)
            % thus the following code converts the measured standard deviation values to a chi square distribution,
            % these values are further transformed to p values (chi2cdf), and then normal values (norminv)
            n = size(sta,4);
            rf = norminv(chi2cdf( n*reshape(rf,[],1).^2 ,n-1));
            
            % set Inf to a reasonable value
            rf(rf==Inf) = 5.5;
            
            % reshape
            rf = reshape(rf,size(sta,1),size(sta,2),size(sta,3));
    end
else
    
    rf = sta;
end

% imagesc(norm_image(rf))
% pause


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZE VALUES
%
% get mean 0, std 1 in each color channel

% reshape to have each stixel's colors in a row
rf_reshape = reshape(rf,[],size(rf,3));

% recenter each color's distribution at 0
%%FIXME
rf_reshape = rf_reshape - repmat( median(rf_reshape) ,size(rf_reshape,1),1);

% get noise sigma of each color channel
noise_sigmas = robust_std(rf_reshape, params.robust_std_method);

% divide by noise sigmas
rf_reshape = rf_reshape./repmat(noise_sigmas,size(rf_reshape,1),1);

% reshape
rf_renorm = reshape(rf_reshape,size(rf,1),size(rf,2),size(rf,3));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE RF STRENGTH
%
% also define function that converts threshold to relevant units
%
% rf_test contains values which will be tested to identify significant stixels
% rf_strength contains values that are output
% f_thresh is the function applied to params.thresh to determine the threshold


% determine how to compute strength
if ischar(params.strength)
    strength_type = params.strength;
elseif iscell(params.strength)
    strength_type = params.strength{1};
end

switch strength_type
    case 'vector length'
        
        % get vector length of each pixel
        rf_strength = sqrt(sum(rf_renorm.^2,3));
        % then the SQUARE of noise pixels will follow a chi square distribution (d.f. = # color dimensions of the original rf)
        
        % rf test is the same as rf strength
        rf_test = rf_strength;
        
        % base threshold on chi square distribution
        f_thresh = @(th)sqrt(chi2inv(normcdf(th),size(rf,3)));
        
    case 'inner'
        
        % make convenient variable for vector to take dot product with 
        v = params.strength{2};
        
        % ensure length of v = number of colors
        if length(v) ~= num_col
            error('Vector for dot product has length (%d) not equal to the number of colors (%d).',length(v),num_col)
        end
        
        % apply dot product to each stixel
        
        % reshape so rf has one color per column
        rf_temp = reshape(rf_renorm,[],num_col);
        
        % take dot product
        rf_dot_prod = rf_temp*v';
            
        % normalize to make sigma 1
        % if the different color channels are independent, this is unnecessary
        % but if they are correlated, then sigma will not be 1
        rf_dot_prod = rf_dot_prod / robust_std(rf_dot_prod, params.robust_std_method);
        
        % reshape
        rf_strength = reshape(rf_dot_prod,size(rf,1),size(rf,2));
        
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
        rf_temp = reshape(rf_renorm,[],num_col);
        
        % go through each vector
        for vv =1:size(vectors,1)
            
            % compute dot product
            temp = rf_temp*vectors(vv,:)';
            
            % normalize to make sigma 1
            % if the different color channels are independent, this is unnecessary
            % but if they are correlated, then sigma will not be 1
            temp = temp/robust_std(temp, params.robust_std_method);
            
            % take dot product AND absolute value
            rf_dot_prod_temp(:,:,vv) = abs(reshape(temp,size(rf,1),size(rf,2)));
            
            % note whether absolute value exceeds threshold
            rf_tests_temp(:,:,vv) = rf_dot_prod_temp(:,:,vv) > -norminv(normcdf(-params.thresh)/2.2/2);
            
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
switch params.select 
    case {'thresh','max'}
        abs_thresh = f_thresh(params.thresh);
end


% look for significant stixels
switch params.select

    case 'thresh' % apply simple threshold

        % return significant stixels
        sig_stixels = rf_test > abs_thresh;
        
        % return rf strength as vector length
        rf_strength_out = rf_strength;


    case 'max'  % return any pixel that is a maximum in its local neighborhood
        
        sig_stixels = find_local_maxima(rf_test,'thresh',abs_thresh,'radius',params.radius,'return','matrix');
        
        % return rf strength as vector length
        rf_strength_out = rf_strength;
        
        
end


% ensure sparse!
sig_stixels = sparse(sig_stixels);

