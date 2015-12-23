function gauss = make_gaussian(varargin)
% MAKE_GAUSSIAN     Compute difference of gaussian values in a matrix or along a line
%
% usage:  gauss = make_gaussian(varargin)
%
% arguments:   params - struct or list of optional parameters (see below)
%
% outputs:      gauss - matrix containing 2D difference of gaussian
%
%
% optional fields in params, their default values, and what they specify:
%
%       center_radius      	1           center gaussian radius
%       center_scale     	1           center gaussian scale
%       surround_radius    	2.5         surround gaussian radius
%       surround_scale      0           surround gaussian scale
%       dim                 2           dimensions in the domain
%       sparse                          return sparse matrix?
%                           true        if effective_radius < Inf
%                           false       if effective_radius = Inf
%
%       normalize           'none'      normalize center and surround gaussian before summing
%                                           'none' - none (center is set to 1)
%                                           'sum'  - set sum to 1
%                                           'var'  - set var to 1
%                                          'norm'  - set L2 norm to 1
%               NOTE: output = cent_scale * norm(center_gaussian) - surr_scale * norm(surr_gaussian)
%
%   NOTE: all parameters are in ABSOLUTE units, i.e. the surround values are NOT relative to the center
%
%  if dim == 2
%       x_size          3           size of domain in x
%       y_size          3           in y
%       center          [2 2]       location of center, domain will be 1:x_size by 1:y_size
%       effective_radius
%                       Inf         effective radius.  defines region in which to compute gaussian.
%                                       outside this region is all 0s.  if Inf, gaussian is computed everywhere.
%
%  if dim == 1
%       x_vals          -5:0.1:5    list of all points in the domain
%       center          0           center point
%   
%
% gauthier 2008-10
%
%
 


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('center_radius', 1);
p.addParamValue('center_scale', 1);
p.addParamValue('surround_radius', 2.5);
p.addParamValue('surround_scale', 0);
p.addParamValue('normalize', 'none', @(x)any(strcmpi(x,{'none','sum','var','norm'})));
p.addParamValue('dim', 2, @(x)any(find(x == [1 2])));

% forked parameters

% add defaults based on number of dimensions
p.KeepUnmatched = true;
if isprop(p, 'PartialMatching'), p.PartialMatching = false; end
p.parse(varargin{:});
switch p.Results.dim

    case 1 % 1D gaussian
        p.addParamValue('x_vals',-5:0.1:5);
        p.addParamValue('center',0,@(x)numel(x) == 1);
    
    case 2 % 2D gaussian
        p.addParamValue('x_size',3);
        p.addParamValue('y_size',3);
        p.addParamValue('center',[2 2],@(x)numel(x) == 2);
        p.addParamValue('effective_radius', Inf);
        
        % add defaults based on effective radius
        p.parse(varargin{:});
        
        
        switch p.Results.effective_radius
            case Inf % computed everywhere, don't bother making sparse
                p.addParamValue('sparse',false);

            otherwise % computed in a small region, make sparse
                p.addParamValue('sparse',true);
        end
end



% resolve user input and default values
p.KeepUnmatched = false;
p.parse(varargin{:});

% get params struct
params = p.Results;



% generate profile function for center and surround

% profile for a single gaussian (used for both center and surround)
% this will be applied to the whole domain, where each entry is the sigma-normalized distance from the center
% thus one_gauss_fcn(1) is the value at 1 sigma
one_gauss_fcn = @(x)( exp(-0.5*(x.^2)) );

% normalization function
% applied to the each gaussian (center and surround) before they are weighted and summed
switch params.normalize
    case 'none'
        norm_gauss = @(g)(g);
    case 'sum'
        norm_gauss = @(g)(g/sum(reshape(g,[],1)));
    case 'var'
        norm_gauss = @(g)(g/var(reshape(g,[],1)));
    case 'norm'
        norm_gauss = @(g)(g/norm(reshape(g,[],1)));
    otherwise
        error('type of normalization ''%s'' not recognized',params.normalize)
end

% combian gaussian profile(s) and normalization functions
%
% if no surround, don't waste computation time by including it
if params.surround_scale == 0
    gauss_fcn = @(x) params.center_scale*norm_gauss( one_gauss_fcn(x/params.center_radius) );
else
    gauss_fcn = @(x) params.center_scale*norm_gauss( one_gauss_fcn(x/params.center_radius) ) - ...
        params.surround_scale*norm_gauss( one_gauss_fcn(x/params.surround_radius) );
end


% % if no surround, don't waste computation time by including it
% if params.surround_scale == 0
%     gauss_fcn = @(x)( params.center_scale*exp(-0.5*(x/params.center_radius).^2) );
% else
%     gauss_fcn = @(x)( params.center_scale*exp(-0.5*(x/params.center_radius).^2) - ...
%         params.surround_scale * exp(-0.5*(x/params.surround_radius).^2) );
% end



% apply profile function to domain

switch params.dim
    case 1
        
        % apply function to x_vals
        gauss = gauss_fcn(params.x_vals - params.center);
        
        
    case 2

        % compute effective region (where to compute the gaussian)
        eff_x = max(1,round(params.center(1))-params.effective_radius) : ...
            min(params.x_size,round(params.center(1))+params.effective_radius);
        eff_y = max(1,round(params.center(2))-params.effective_radius) : ...
            min(params.y_size,round(params.center(2))+params.effective_radius);
        eff_size_x = length(eff_x);
        eff_size_y = length(eff_y);
        % re-compute center in effective region
        eff_ctr = params.center - [min(eff_x) min(eff_y)] + 1;


        % compute distance of each point from the center
        x_dist = ((1:eff_size_x) - eff_ctr(1)).^2;
        y_dist = ((1:eff_size_y) - eff_ctr(2)).^2;
        dist = sqrt(ones(eff_size_y,1)*x_dist + y_dist'*ones(1,eff_size_x));

        % apply profile to these distances
        gauss_eff = gauss_fcn(dist);

        % make matrix of zeros
        if params.sparse
            gauss = sparse(params.y_size,params.x_size);
        else
            gauss = zeros(params.y_size,params.x_size);
        end
            
        % put gaussian shape in the effective region
        gauss(eff_y,eff_x) = gauss_eff;
end

