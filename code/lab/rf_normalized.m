function  [rf_norm,extras,params] = rf_normalized(rf, varargin)
% rf_normalized     Normalize a RF based on peak, variance in a ROI, etc
%
% usage:  [rf_norm,extras,norm_params] = rf_normalized(rf, params)
%
% arguments:       rf - STA spatial rf, 2-d or 3-d matrix (height, width, color)
%            varargin - struct or list of optional parameters (see below)
%
%
% outputs:         rf_norm - normralized RF
%                   extras - struct of additional information
%                   params - copy of the parameters whih were used
%
%
% optional params, their default values, and what they specify:
%
% method            'variance'  type of normalization to perform
%                                	'variance' - set each RF to have unit variance within a ROI
%                                   'peak'     - set the largest stixel to a value of 1
%
%                         
%  if method == 'variance'
%       roi         true([size(rf,1) size(rf,2)]))
%                               boolean matrix of which stixels to include in the region of interest
%                                   if not specified, all stixels are used, which is ridiculous
%                         
%
%  if method == 'peak'
%       peak      	1         	value of the peak in the normalized RF
%
% 2009-04  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('method','variance', @(x)any(strcmpi(x,{'variance','peak'})));

% add other parameters based on 'method'
p.KeepUnmatched = true;
p.parse(varargin{:});
switch p.Results.method
    case 'variance'
        p.addParamValue('roi', true([size(rf,1) size(rf,2)]));
    case 'peak'
        p.addParamValue('peak', 1);
end
p.KeepUnmatched = false;

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;





switch params.method
    
    case 'variance' % normalize variance within the ROI
        
        % get values in the ROI
        
        % if they're all included, just get them all
        if all(all(params.roi == true([size(rf,1) size(rf,2)])))
            roi_vals = reshape(rf,[],1);
            
        else
            % get values at these indices, from each color
            roi_vals = [];
            for cc = 1:size(rf,3)
                tmp = rf(:,:,cc);
                roi_vals = [roi_vals reshape(tmp(params.roi),[],1)];
            end
        end
        
        % divide the rf by the std of values in the ROI
        rf_norm = rf/std(roi_vals);
        
        % no extra info
        extras = struct;
    
    case 'peak'
        % 0 -> 0, peak -> 1
        
        % get peak
        peak = max(rf(:));
        
        % divide by it
        rf_norm = rf / peak;
        
        % note the peak 
        extras.peak = peak;
        
    otherwise
        error('rf_normalized: type of normalization not recognized (%s).',params.norm_type)
        
end
