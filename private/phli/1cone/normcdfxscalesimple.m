function [p resnorm residual h] = normcdfxscalesimple(crs, crsx, varargin)
% NORMCDFXSCALE
% usage: p = normcdfxscale(crs, crsx, varargin);
%
% Simplified version of NORMCDFXSCALE with only negative contrasts 
% (primarily for allcones analysis)
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

allmapindices = 1:size(crsx,1);

opts = inputParser;
opts.addParamValue('subtract_baseline', false);
opts.addParamValue('baseline', []);
opts.addParamValue('fitfunc', @normcdffitn);
opts.addParamValue('MaxFunEvals', 5000);
opts.addParamValue('MaxIter', 5000);
opts.addParamValue('optimset', {});
opts.addParamValue('plot', false);
opts.addParamValue('Markers', []);
opts.addParamValue('MarkerSize', 15);
opts.addParamValue('fill', false);
opts.addParamValue('colors', []);
opts.addParamValue('cmf', @jet);
opts.addParamValue('p0', []);
opts.addParamValue('sat0', []);
opts.addParamValue('title', true);
opts.parse(varargin{:});
opts = opts.Results;

optimopts = optimset('MaxFunEvals', opts.MaxFunEvals, 'MaxIter', opts.MaxIter, opts.optimset{:});
n = size(crsx,1);

% Subtract baseline?
if opts.subtract_baseline
    
    % Assumes baseline is just location where CRSX == 0.  There may be one 
    % of these for each row, but assumes that it is actually the same value 
    % copied to each row.
    if isempty(opts.baseline)
        blanks = crsx(:) == 0;
        baselines = crs(blanks);
        opts.baseline = baselines(1);
    end
    
    crs = crs - opts.baseline;
end


% Guess reasonable initial conditions if not given
if isempty(opts.p0)
    if isempty(opts.sat0)
        opts.sat0 = max(crs(:));
    end
    sigma0 = 0.5;
    opts.p0 = [zeros(1,size(crs,1)) opts.sat0 sigma0];
end

% Set fitting bounds
lb = zeros(1,n+2);
ub = Inf(1,n+2);

% Run fit
[p resnorm residual] = lsqcurvefit(opts.fitfunc, opts.p0, crsx, crs, lb, ub, optimopts);

 
% Plot?
if opts.plot
    oldhold = ishold();
    if oldhold, cla; end
    hold on;

    % Undo baseline subtraction
    if opts.subtract_baseline
        crs = crs + opts.baseline;
    end

    if isempty(opts.colors)
        opts.colors = opts.cmf(n);
    end
        
    for i = 1:n
        x = crsx(i,:) .* p(i);
%         h(i) = plot(x, crs(i,:), '.', 'Color', opts.colors, 'MarkerSize', opts.MarkerSize);
        h(i) = plot(x, crs(i,:), '.', 'Color', opts.colors(i,:), 'MarkerSize', opts.MarkerSize);
    end
    
    if ~isempty(opts.Markers)
        for i = 1:n
            set(h(i), 'Marker', opts.Markers(i));
        end
    end
    
    if opts.fill
        for i = 1:n
            set(h(i), 'MarkerFaceColor', opts.colors(i,:));
        end
    end
    
    if opts.title
        title(['Params: ' num2str(p) '  Error: ' num2str(resnorm)]);
    end
    
    if ~oldhold, hold off; end
end


if nargout == 0, clear p; end


function y = normcdffitn(p,x)
n = size(x,1);
sat   = p(n+1);
sigma = p(n+2);
for i = 1:n
    xscale = p(i);
    y(i,:) = sat .* normcdf(abs(x(i,:)) * xscale, 1, sigma);
end