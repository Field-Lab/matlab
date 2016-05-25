function [threshold] = sig_stixels_threshold(sta, false_stixels, varargin)
% sig_stixels _threshold    identify threshold to input into
% significant_stixels function
%
% strategy:
%
%  * spatially filter the STA (if desired)
%  * collapse across time
%  * reduce to 1 color
%  * compute where the threshold should be on the null distribution to
%  achieve the inputted false positive stixels on the true distribution
%
%
%
% usage:  [threshold] = sig_stixels_threshold(sta, false_stixels, params)
% Then call sig_stixels =significant_stixels(sta, 'select', 'thresh', 'thresh', threshold)
% The sta inputted to sig_stixels_threshold is computed from the null
% distribution and the sta inputted to significant_stixels is computed from
% the true movie file
%
% arguments:     sta - 2d, 3d, or 4d matrix (height, width, color, time)
%                   computed from a movie file with the wrong seed
%                false_stixels - the number of stixels that contain noise
%                   that would be chosen as significant stixels
%                   i.e. if you choose 0.5 false stixels, every other time
%                   significant stixels are calculated, 1 stixel in the
%                   noise will be selected as significant
%                params - struct of optional parameters (see below)
%
% outputs:       threshold - the value to input into sig_stixels
%                   function in order to calculate sig_stixels with the selected false
%                   positive number of stixels in each STA
%
%
% optional fields in params, their default values, and what they specify:
%
% filter      	[]          filter to convolve with each color channel in each frame
%
% robust_std_method 3       which algorithm to use for calculating
%                               robust_std; see robust_std for more
%                               information.
%
% 2015-04-27 Colleen Rhoades
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('filter', []);
p.addParamValue('robust_std_method', 3);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% if STA is empty, return not a number
if isempty(sta)
    threshold = Nan;
    return
end

%% APPLY SPATIAL FILTER

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

% reshape to have each stixel's colors in a row
sta_reshape = reshape(permute(sta_filt,[1 2 4 3]),[],size(sta,3));

% get noise sigma of each color channel
noise_sigmas = robust_std(sta_reshape, params.robust_std_method);

% divide by noise sigmas
sta_reshape = sta_reshape./repmat(noise_sigmas,size(sta_reshape,1),1);

% reshape
sta_renorm = permute(reshape(sta_reshape,size(sta,1),size(sta,2),size(sta,4),size(sta,3)),[1 2 4 3]);

rf = std(sta_renorm,1,4);

%% transform to normal distribution

% the distribution of standard deviation values with 1/n normalization (s_hat) has the following property:
%   (n-1*)(s_hat^2)/(s_true^2)  =  chi square with n-1 degrees of freedom
%       where n is the number of samples and s_true is the true standard deviation
% here s_true is 1 (by the normalization above), and n = size(sta,4)
% thus the following code converts the measured standard deviation values to a chi square distribution,
% these values are further transformed to p values (chi2cdf), and then normal values (norminv)
n = size(sta,4);
rf = norminv(chi2cdf((n-1)*reshape(rf,[],1).^2 ,n-1));

rf(isinf(double(rf))) = 5.5;
% Fit a normal distribution to the distribution of stixel values
Pd = fitdist(double(rf),'normal'); 

% 1/2 the false stixels in each tail
p = false_stixels/2/length(rf); 
threshold = Pd.mean - norminv(p, Pd.mean, Pd.sigma);

