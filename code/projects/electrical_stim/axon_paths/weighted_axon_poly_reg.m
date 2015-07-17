function [curve_x, curve_y, p, soma_x, soma_y, valid, res] = weighted_axon_poly_reg(eiAmps, varargin)
% AXON_POLY_REG() takes in an EI and finds a polynomial
% line of best fit for the weighted amplitudes.
%   inputs:        eiAmps
%           end_electrode - The endpoint of the axon approximation. Soma
%                           location is estimated, but end electrod must
%                           be specified by number. (ex: 509 out of 512)
% 
%   optional:        plot - If 'plot' is followed by 'true',the
%                              function plots the electrodes and curve.
%                       N - Specifies the degree of polynomial. 
%               ei_thresh - Manually set the ei threshold of points
%                           considered for the axon curve.
%   outputs:      curve_x - X coordinates of axon estimation.
%                 curve_y = Y coordinates of axon estimation.

% Yields the minimum squared-error, on average
N = 7;


% Thresholds that yield the best fit vary from cell to cell.
% 5 seems to be a good comproimse.
point_threshold = max(eiAmps) / 20;

plot_reg = false;
valid = true;

% Increasing this would seem to give more fine differentiation between
% ei amps, but increasing it doesn't decrease error overall. 
amp_scaling = 1;

% Read if plot is true.
nbin = length(varargin);
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch varargin{j*2-1}
        case 'plot'
            plot_reg = varargin{j*2};   
        case 'N'
            N = varargin{j*2};
        case 'ei_thresh'
            point_threshold = varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% Get electrode coordinates
switch length(eiAmps)
    case 512
        [coords(:,1),coords(:,2)] = getElectrodeCoords512();
    case 519
        [coords(:,1),coords(:,2)] = getElectrodeCoords519();
    case {61,64}
        [coords(:,1),coords(:,2)] = getElectrodeCoords61();
    otherwise
        err = MException('MATLAB:InvArgIn',...
            'Unknown array specified - must be 512, 519, or 61');
        throw(err);
end


% Finds array bounds for creating a grid and plotting.
arrayX = max(abs(coords(:,1)));
arrayY = max(abs(coords(:,2)));


% Find soma location as in eiContour_wLinFit
linFitThresh = 6;
[~,c] = size(eiAmps);
if c == 1
    eiAmps = eiAmps';
end   
[~,col] = find(eiAmps > linFitThresh);
aa = round(eiAmps(col))';
yy = coords(col,2)';
xx = coords(col,1)';
sortaa = sort(aa,1,'descend')';
largestAmps = sortaa(1:2);
[~,IA,~] = intersect(aa,largestAmps);

aa = aa';
soma_x = 1/sum(largestAmps) * sum(xx(IA).*aa(IA));
soma_y = 1/sum(largestAmps) * sum(yy(IA).*aa(IA));

% Sets the conversion from ei amplitude to repetition of xy coordinate
counts = floor(eiAmps.*amp_scaling);

above_thresh_points = [];
above_thresh_amps = [];

% The regression will only consider points above the point threshold, which
% are plotted repeatedly in proportion to their amplitude.

% Locations of the electrodes above point_threshold
above_thresh_x = [];
above_thresh_y = [];

for n = 1:length(eiAmps)
    point = coords(n, :);
    if counts(n) > (point_threshold * amp_scaling)
        above_thresh_points = vertcat(above_thresh_points, point);
        above_thresh_amps = vertcat(above_thresh_amps, eiAmps(n));
        above_thresh_x = vertcat(above_thresh_x, point(1));
        above_thresh_y = vertcat(above_thresh_y, point(2));
    end    
end

if size(above_thresh_x) < 20
    valid = false;
    warning('weighted_axon_poly_reg:Unclear',...
    'Threshold is too high, or there are too few points with signal for a good fit. ');
end

% If the axon path seems like it might be not a function (if it may be
% aligned vertically- y spread greater than x), then x and y will be
% swapped prior to finding the line of best fit. 
if range(above_thresh_y) > range(above_thresh_x)
    swap = above_thresh_points(:,2);
    above_thresh_points(:,2) = above_thresh_points(:,1);
    above_thresh_points(:,1) = swap;
    original_soma_x = soma_x;
    soma_x = soma_y;
    first_coord = min(above_thresh_y);
    last_coord = max(above_thresh_y);
    search_coords = coords(:,2);
    swap = arrayY;
    arrayY = arrayX;
    arrayX = swap;
else
    first_coord = min(above_thresh_x);
    last_coord = max(above_thresh_x);
    search_coords = coords(:,1);
end    
    

% Creates a polynomial line of degree N that best fits the EI amps 
% with repetitions to represent strength of signal
%degree = strcat('poly',num2str(N));
%f = fittype(degree);
%options=fitoptions(degree);
%options.Weights = eiAmps;
%[fun,p]=fit(coords(:,1),coords(:,2),f,options);

% Can do the same thing without checking thesholds, because the polynomial
% is weighted, the below-threshold points barely affect the fit

% x = coords(:,1);
% y = coords(:,2);
% p = polyfitweighted(x, y, N, eiAmps);

x = above_thresh_points(:,1);
y = above_thresh_points(:,2);


% Decides which way the axon should go, based on the average value of the
% electrodes to the right and left of the soma. 

left_of_soma = find(search_coords < soma_x);

indicies = 1:length(eiAmps);
right_of_soma = setxor(indicies, left_of_soma);

% Pad to prevent a small range from having artificially high average
% values.
val_left = sum(eiAmps(left_of_soma));
val_right = sum(eiAmps(right_of_soma));

x_steps = 1;

% Sets range for axon x
if val_left > val_right
    curve_x = first_coord:x_steps:soma_x;
    curve_x = fliplr(curve_x);
    x_range = soma_x - first_coord;
else
    curve_x = soma_x:x_steps:last_coord;
    x_range = last_coord - soma_x;
end  

soma_range = floor((arrayX*2)/x_range);

%if soma_range > 7
%    valid = false;
    %warning('weighted_axon_poly_reg:Short',...
    %    'Threshold is too high, or axon path is too small for a valid fit');
%else    
    while soma_range > 1
        N = N - 1;
        soma_range = soma_range - 1;
    end    
%end 

if N < 3
    N = 3;
end    

p = polyfitweighted(x, y ,N,above_thresh_amps);
curve_y = polyval(p,curve_x);

if abs(max(curve_y)) > (arrayY + 50)
    valid = false;
    warning('weighted_axon_poly_reg:OutOfBounds',...
    'Axon estimation goes out of array bounds. This is usually the result of an EI with too few electrodes with signal.' );
end    

k = find(abs(curve_y) > arrayY);
if ~isempty(k)
    curve_x = curve_x(1:min(k));
    curve_y = curve_y(1:min(k));
end
% If the x and y were switched before the regression, switch them back.
if range(above_thresh_y) > range(above_thresh_x)
    swap = curve_x;
    curve_x = curve_y;
    curve_y = swap;
    soma_x = original_soma_x;
end    

res = 0;
for elec = 1:length(eiAmps)
    min_dist_to_axon = sqrt( min( sum( bsxfun(@minus, vertcat(curve_x, curve_y), [coords(elec,1);coords(elec,2)]).^2,1)));
    res = res + (min_dist_to_axon * eiAmps(elec));
end
res = res/max(eiAmps);

if res > 3000
    valid = false;
    warning('weighted_axon_poly_reg:NoDefinedAxon', ...
        'EI is too scattered for a clear axon path');


if plot_reg
    figure
    scatter(coords(:,1),coords(:,2), eiAmps+.1, 'filled');
    hold on;
    plot(curve_x, curve_y, '-');
end   

if ~valid
    curve_x = [];
    curve_y = [];
end    

end

