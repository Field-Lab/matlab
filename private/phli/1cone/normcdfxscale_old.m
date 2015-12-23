function [p resnorm residual] = normcdfxscale_old(crs, crsx, varargin)
% NORMCDFXSCALE
% usage: p = normcdfxscale(crs, crsx, varargin);
%
% In this version, for some reason I thought we could use a single
% saturation value for both negative and positive CRSX.  Not sure why I
% thought that would make sense.  Maybe a little more for a situation where
% we're calculating the C/R with different pos/neg templates, but still
% doesn't really make sense.
%
% CRS is an MxN matrix where M is the number of cones and N is the number
% of contrasts.  CRSX is also an MxN matrix.
%
% Because of the way lsqcurvefit works, each row must have the same number
% of data points.  If the data do not have the same x locations for each
% row, the input matrices have to be padded with NaNs (?).
%
% 2012-02 phli
%

opts = inputParser;
opts.addParamValue('rectfit', true);
opts.addParamValue('polarity', 1);
opts.addParamValue('fitfunc', []);
opts.addParamValue('MaxFunEvals', 5000);
opts.addParamValue('MaxIter', 5000);
opts.addParamValue('optimset', {});
opts.addParamValue('plot', false);
opts.addParamValue('MarkerSize', 15);
opts.addParamValue('colors', []);
opts.addParamValue('cmf', @hsv);
opts.addParamValue('p0', []);
opts.addParamValue('sat0', []);
opts.addParamValue('negxscale0', 0.5);
opts.parse(varargin{:});
opts = opts.Results;


if isempty(opts.fitfunc)
    if opts.rectfit, opts.fitfunc = @normcdfrectfitn;
    else opts.fitfunc = @normcdffitn; end
end


optimopts = optimset('MaxFunEvals', opts.MaxFunEvals, 'MaxIter', opts.MaxIter, opts.optimset{:});
n = size(crsx,1);


if isempty(opts.sat0)
    opts.sat0 = max(abs(crs(:)));
end
if isempty(opts.p0)
    sat0 = opts.sat0;
    sigma0 = 0.5;
    opts.p0 = [ones(1,n) sat0 sigma0 opts.negxscale0];
end


lb = [zeros(1,n) -Inf 0 0]; % Let sat be anything to accomodate off cells
ub = [];
[p resnorm residual] = lsqcurvefit(opts.fitfunc, opts.p0, crsx, crs, lb, ub, optimopts);


if opts.plot
    oldhold = ishold();
    if oldhold, cla; end
    hold on;
    
    if isempty(opts.colors)
        opts.colors = opts.cmf(n);
    end
    
    for i = 1:n
        plot(crsx(i,:).*p(i), crs(i,:), '.', 'Color', opts.colors(i,:), 'MarkerSize', opts.MarkerSize);
    end
    title(['Params: ' num2str(p) '  Error: ' num2str(resnorm)]);
    
    if ~oldhold, hold off; end
end

if nargout == 0, clear p; end


% Same as below but rectified
function y = normcdfrectfitn(p,x)
n = size(x,1);

sat       = p(n+1);
sigma     = p(n+2);
negxscale = p(n+3); % Relative to positive

for i = 1:n
    xscale = p(i);
    
    xpos = x(i,:) >= 0;
    xneg = ~xpos;
    
    y(i,xpos) = sat .* normcdf(               x(i,xpos) * xscale, 1, sigma);
    y(i,xneg) = sat .* normcdf(negxscale.*abs(x(i,xneg))* xscale, 1, sigma);
end


% Gaussian CDF fit
function y = normcdffitn(p,x)
n = size(x,1);

sat       = p(n+1);
sigma     = p(n+2);
negxscale = p(n+3); % Relative to positive

for i = 1:n
    xscale = p(i);
    
    xpos = x(i,:) >= 0;
    xneg = ~xpos;
    
    y(i,xpos) =  sat .* normcdf(               x(i,xpos)  * xscale, 1, sigma);
    y(i,xneg) = -sat .* normcdf(abs(negxscale.*x(i,xneg)) * xscale, 1, sigma);
end