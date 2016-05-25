function [curve_x, curve_y, p, soma_x, soma_y, valid, res] = weighted_axon_poly_reg(eiAmps, varargin)
% AXON_POLY_REG() takes in an EI and finds a polynomial
% line of best fit for the weighted amplitudes.
%   inputs:        eiAmps - either the max ei amps (512 x 1) or the entire
%                           EI may be imputted. The entire EI enables
%                           for possible axon speed/temporal testing, but
%                           needs more work.
%           end_electrode - The endpoint of the axon approximation. Soma
%                           location is estimated, but end electrod must
%                           be specified by number. (ex: 509 out of 512)
% 
%   optional:        plot - If 'plot' is followed by 'true',the
%                              function plots the electrodes and curve.
%                       N - Specifies the degree of polynomial. 
%               ei_thresh - Manually set the ei threshold of points
%                           considered for the axon curve.
%           axonBundleFit - default false, set to true to ignore the
%                           'soma' and fit a across the entire array
%             newXYCoords - alternate set of coordinates to use for
%                           plotting fits over an electrode array.
%   outputs:      curve_x - X coordinates of axon estimation.
%                 curve_y - Y coordinates of axon estimation.
%                       p - polynomial fit coefficients
%                   valid - boolean of whether or not the fit is valid
%                     res - residual, measures goodness of fit
% Alena Rott, summer 2015
% LG added axon bundle fitting, ability to plot with arbitrary coordinates.

% Yields the minimum squared-error, on average
N = 7;
ei = [];

% Thresholds that yield the best fit vary between cells and datasets
point_threshold = max(eiAmps) / 15;

% Set up default parameters.
plot_reg = false;
valid = true;
axonBundleFit = false; 
newXYCoords = []; 

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
        case 'axonBundleFit'
            axonBundleFit = varargin{j*2};
        case 'newXYCoords'
            newXYCoords = varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

hexagonal = false;

% Get electrode coordinates
switch length(eiAmps)
    case 512
        [coords(:,1),coords(:,2)] = getElectrodeCoords512();
        load('adj_mat_512.mat');
        adj_matrix = adj_mat_512;
    case 519
        [coords(:,1),coords(:,2)] = getElectrodeCoords519();
        load('adj_mat_519.mat');
        adj_matrix = adj_mat_519;
        hexagonal = true;
    case {61,64}
        [coords(:,1),coords(:,2)] = getElectrodeCoords61();
        load('adj_mat_61.mat');
        adj_matrix = adj_mat_61;
    otherwise
        err = MException('MATLAB:InvArgIn',...
            'Unknown array specified - must be 512, 519, or 61');
        throw(err);
end

if ~isempty(newXYCoords)
    coords = newXYCoords;
end
% standardizes eiAmps to a column vector
if size(eiAmps, 1) == 1
    eiAmps = eiAmps';
end 

if size(eiAmps, 2) ~= 1
    ei = eiAmps;
    eiAmps = max(ei) - min(ei);
    if size(eiAmps, 1) == 1
        eiAmps = eiAmps';
    end 
    
    if size(ei, 1) == 71
        ei = ei';
    end;    
end    

% Finds array bounds for creating a grid and plotting.
arrayX = max(abs(coords(:,1)));
arrayY = max(abs(coords(:,2)));

[max_amp,max_amp_i] = max(eiAmps);
adj = adj_matrix{max_amp_i};
[max_adj,max_adj_i] = max(eiAmps(adj));

indicies = [max_amp_i;adj(max_adj_i)];

% Soma location is center of mass of largest amplitude and largest neighbor
soma_x = sum(coords(indicies,1).*eiAmps(indicies)) / (max_amp + max_adj);
soma_y = sum(coords(indicies,2).*eiAmps(indicies)) / (max_amp + max_adj);



above_thresh_points = [];
above_thresh_amps = [];

% Locations of the electrodes above point_threshold
above_thresh_x = [];
above_thresh_y = [];

for n = 1:length(eiAmps)
    point = coords(n, :);
    if eiAmps(n) > point_threshold
        above_thresh_points = vertcat(above_thresh_points, point);
        above_thresh_amps = vertcat(above_thresh_amps, eiAmps(n));
        above_thresh_x = vertcat(above_thresh_x, point(1));
        above_thresh_y = vertcat(above_thresh_y, point(2));
    end    
end

if size(above_thresh_x) < 18
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
    switched = true;
else
    first_coord = min(above_thresh_x);
    last_coord = max(above_thresh_x);
    search_coords = coords(:,1);
    switched = false;
end    

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

% Signal on either side of soma for deciding which direction the axon goes.
val_left = sum(eiAmps(left_of_soma))*abs(soma_x - first_coord);
val_right = sum(eiAmps(right_of_soma))*abs(soma_x - last_coord);

x_steps = 1;

if axonBundleFit
    % Sets range for axon x
    if val_left > val_right
        curve_x = first_coord:x_steps:last_coord;
        curve_x = fliplr(curve_x);
        x_range = last_coord - first_coord;
    else
        curve_x = first_coord:x_steps:last_coord;
        x_range = last_coord - first_coord;
    end
else
    % Sets range for axon x
    if val_left > val_right
        curve_x = first_coord:x_steps:soma_x;
        curve_x = fliplr(curve_x);
        x_range = soma_x - first_coord;
    else
        curve_x = soma_x:x_steps:last_coord;
        x_range = last_coord - soma_x;
    end
end
% The following code that alters N is meant to give smaller axons an
% estimate with a smaller-order polynomial (smaller slices of the array
% don't need a 7th order polynomial)
soma_range = floor((arrayX*2)/x_range);
   
while soma_range > 1
    N = N - 1;
    soma_range = soma_range - 1;
end    

% The order should be at least 3
if N < 3
    N = 3;
end    

p = polyfitweighted(x, y ,N,above_thresh_amps);
curve_y = polyval(p,curve_x);

if max(abs(curve_y)) > (arrayY + 50)
    valid = false;
    warning('weighted_axon_poly_reg:OutOfBounds',...
    'Axon estimation goes out of array bounds. This is usually the result of an EI with too few electrodes with signal.' );
end    

test = [curve_x; curve_y];

% indicies of where the curve goes off the array
% (a bit complicated for hexagonal array)
k = [];
if valid && ~switched && hexagonal
    for h = 1:length(curve_x)
        if abs(curve_y(h)) + (5/12)*curve_x(h) > 360
            k = [k h];
        end    
    end
elseif valid && hexagonal
    for h = 1:length(curve_x)
        if abs(curve_y(h)) + (12/5)*abs(curve_x(h)) > 864
            k = [k h];
        end    
    end
else    
    k = find(abs(curve_y) > arrayY);
end    
if ~isempty(k)
    %if hexagonal
    %    curve_x = curve_x(max(k):length(curve_x));
    %    curve_y = curve_y(max(k):length(curve_y));
    %else    
        curve_x = curve_x(1:min(k));
        curve_y = curve_y(1:min(k));
    %end    
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

if res > 10000
    valid = false;
    warning('weighted_axon_poly_reg:NoDefinedAxon', ...
        'EI is too scattered for a clear axon path');
end

if plot_reg
    figure
    scatter(coords(:,1),coords(:,2), (eiAmps+0.1)*8, 'filled', 'MarkerFaceColor', 'black');
    hold on;
    plot(curve_x, curve_y, 'r');
    scatter(soma_x, soma_y, max(eiAmps), 'MarkerFaceColor', 'red'); 
end   

% Amacrine testing- when the whole EI is included. Currently doesn't work
if ~isempty(ei)
    spike = [];
    
    for t = 1:size(ei, 2)
        spike = [spike; sum(ei(:,t))];
    end 
    
    [~,spike_max] = max(spike);
    if spike_max > 35
        warning('weighted_axon_poly_reg:PossibleAmacrine');
        valid = false;
    end
end


%if ~valid
%    curve_x = [];
%    curve_y = [];
%end 


end

