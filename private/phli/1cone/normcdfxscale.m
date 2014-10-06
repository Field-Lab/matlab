function [p resnorm residual] = normcdfxscale(crs, crsx, varargin)
% NORMCDFXSCALE
% usage: p = normcdfxscale(crs, crsx, varargin);
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
opts.addParamValue('subtract_baseline', false);
opts.addParamValue('baseline', []);
opts.addParamValue('fitfunc', @normcdffitn);
opts.addParamValue('MaxFunEvals', 5000);
opts.addParamValue('MaxIter', 5000);
opts.addParamValue('optimset', {});
opts.addParamValue('plot', false);
opts.addParamValue('MarkerSize', 15);
opts.addParamValue('colors', []);
opts.addParamValue('cmf', @jet);
opts.addParamValue('p0', []);
opts.addParamValue('psat0', []);
opts.addParamValue('nsat0', []);
opts.addParamValue('negxscale0', 1);
opts.addParamValue('title', true);
opts.addParamValue('plotfit', true);
opts.parse(varargin{:});
opts = opts.Results;

optimopts = optimset('MaxFunEvals', opts.MaxFunEvals, 'MaxIter', opts.MaxIter, opts.optimset{:});
n = size(crsx,1);
polarity = opts.polarity;


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
    if isempty(opts.psat0)
        xpos = crsx(:) >= 0;
        
        if opts.rectfit
            opts.psat0 = max(crs(xpos));
        else
            % Use polarity to correct
            opts.psat0 = polarity .* max(crs(xpos).*polarity);
        end
    end
    
    if isempty(opts.nsat0)
        xneg = crsx(:) <= 0;
        if opts.rectfit
            opts.nsat0 = max(crs(xneg));
        else
            % Use polarity to correct
            opts.nsat0 = polarity .* min(crs(xpos).*polarity);
        end
    end
    
    psat0 = opts.psat0;
    nsat0 = opts.nsat0;
    sigma0 = 0.5;
    opts.p0 = [ones(1,n) psat0 nsat0 sigma0 opts.negxscale0];
end


% Set fitting bounds; depends on the polarity/rectification situation
if opts.rectfit
    % All constrained to positive
    lb = zeros(1,n+4);
    ub = [];
else
    if polarity == 1
        lb = [zeros(1,n)   0 -Inf   0   0];
        ub = [Inf(1,n)   Inf    0 Inf Inf];
    else
        lb = [zeros(1,n) -Inf   0   0   0];
        ub = [Inf(1,n)      0 Inf Inf Inf];
    end
end


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
        xneg = x < 0;
        
        % Plot with different negative x-scale factor?
        % x(xneg) = x(xneg) .* p(end);

        plot(x, crs(i,:), '.', 'Color', opts.colors(i,:), 'MarkerSize', opts.MarkerSize);
    end
    
    if opts.plotfit
        scaledxranges = minmax(crsx) .* [p(1:4)' p(1:4)'];
        scaledxrange = minmax(scaledxranges(:)');
        
        psat      = p(n+1);
        nsat      = p(n+2);
        sigma     = p(n+3);
        negxscale = p(n+4);
        
        negx = scaledxrange(1):0.1:0;
        negy = nsat .* normcdf(abs(negx) .* negxscale, 1, sigma);
        posx = 0:0.1:scaledxrange(2);
        posy = psat .* normcdf(posx, 1, sigma);

        if opts.subtract_baseline, 
            negy = negy + opts.baseline;
            posy = posy + opts.baseline;
        end
        plot(negx, negy);
        plot(posx, posy);
    end
    
    if opts.title
        title(['Params: ' num2str(p) '  Error: ' num2str(resnorm)]);
    end
    
    if ~oldhold, hold off; end
end


if nargout == 0, clear p; end


function y = normcdffitn(p,x)
n = size(x,1);

psat      = p(n+1);
nsat      = p(n+2);
sigma     = p(n+3);
negxscale = p(n+4); % Relative to positive

for i = 1:n
    xscale = p(i);
    
    xpos = x(i,:) >= 0;
    xneg = ~xpos;
    
    y(i,xpos) = psat .* normcdf(               x(i,xpos) * xscale, 1, sigma);
    y(i,xneg) = nsat .* normcdf(negxscale.*abs(x(i,xneg))* xscale, 1, sigma);
end