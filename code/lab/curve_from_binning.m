function [average_x, average_y, error_y, bin_counts, bin_edges] = curve_from_binning(x,y, varargin)
% curve_from_binning     Interpolate a curve through an x-y scatter of points by
%                           binning on the x axis and computing the "average" x and y values in each bin,
%                           where the "average" can be computed in several ways
%                           
%
% usage:  [average_x, average_y, error_y, bin_counts, bin_edges] = curve_from_binning(x,y, varargin)
%
%
% arguments:      x,y - vectors of the data
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    average_x - "average" x value in the bin (see options below)
%             average_y - "average" y value in the bin (see options below)
%               error_y - "error" on the estimate of the y value (see options below)
%            bin_counts - number of points in each bin
%             bin_edges - vector specifying bin edges
%
%
%
%   NOTE: when x and y are cell arrays of the same length, the outputs are also cell arrays
%
%
% optional params, their default values, and what they specify:
%
% bin_edges        	min(x)+[0:.1:1]*(max(x)-min(x))
%                               edges of the bins.  must be in ascending order.
% num_bins          []          integer, how many bins to use
%                                   if empty, bin_edges is used to determine the bins
%                                   if nonempty, bins are created which each contain an equal number of points
% average_y         'median'    how to compute the "average" y value
%                                   'mean'   - 
%                                   'median' - 
%                                   'sum'    - 
% average_x         'mean'      how to compute the "average" x value
%                                   'mean'   - 
%                                   'median' - 
%                                   'bin center' - use the center point of each bin
% error_y           'std_err'   how to compute error bars on the y value
%                                   'std_err' - standard error, i.e. std/sqrt(n)
%                                   'std'     - standard deviation
%
%
% gauthier  2008-10
% gauthier  2009-09  added 'bin center','num_bins' options
%




% check for multiple lists
if iscell(x)
    % for each element of the cell array
    for cc=1:length(x)
        % bin the data
        [average_x{cc}, average_y{cc}, error_y{cc}, bin_counts{cc}] = curve_from_binning(x{cc},y{cc},varargin{:});
    end
    return
end



% if empty, return empties
if isempty(x)
    average_x = [];
    average_y = [];
    error_y = [];
    bin_counts = [];
    return
end



% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('bin_edges', min(x)+[0:.1:1]*(max(x)-min(x)), @(x)(numel(x)==length(x)) && all(x==sort(x)));
p.addParamValue('num_bins', [],@(x)mod(x,1)==0 && x > 0);
p.addParamValue('average_y', 'median',@(x)any(strcmpi(x,{'mean','median','sum'})));
p.addParamValue('average_x', 'mean',@(x)any(strcmpi(x,{'mean','median','bin center'})));
p.addParamValue('error_y', 'std_err',@(x)any(strcmpi(x,{'std_err','std'})));


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;
%keyboard





% BODY OF THE FUNCTION

% ensure x and y are column vectors
if size(x,2)>1; x = x'; end
if size(y,2)>1; y = y'; end

% set bin edges
if isempty(params.num_bins)
    % if bin count is not specified, use params.bin_edges
    bin_edges = params.bin_edges;
else
    % if bin count is specified, use it
    
    % sort x values
    x_sort = sort(x);
    % take evenly spaced points
    bin_edges(1) = x_sort(1);
    for bb=1:params.num_bins
       bin_edges(bb+1) = x_sort(max(1,floor(bb*length(x)/params.num_bins)));
    end
    
end

% note number of bins
num_bins = length(bin_edges)-1;

% get indices of x that fall within each bin
for bb =1:num_bins
    bin_indices{bb} = find( (x>=bin_edges(bb)) & (x<bin_edges(bb+1)) );
end
% put points that are equal to the largest bin edge in the final bin
bin_indices{num_bins} = [bin_indices{num_bins}; find(x == bin_edges(num_bins+1))];


% AVERAGE X VALUES

% set algorithm
switch params.average_x
    case 'mean'
        f_avg_x = @(nums)mean(nums);
    case 'median'
        f_avg_x = @(nums)median(nums);
end

% fill in variable
if ~strcmp(params.average_x,'bin center')
    average_x = zeros(num_bins,1);
    for bb=1:num_bins
        average_x(bb) = f_avg_x(x(bin_indices{bb}));
    end
else
    average_x = zeros(num_bins,1);
    for bb=1:num_bins
        average_x(bb) = mean(bin_edges([bb bb+1]));
    end
end


% AVERAGE Y VALUES
if nargout >=2
    
    % set algorithm
    switch params.average_y
        case 'mean'
            f_avg_y = @(nums)mean(nums);
        case 'median'
            f_avg_y = @(nums)median(nums);
        case 'sum'
            f_avg_y = @(nums)sum(nums);
    end

    % fill in variable
    average_y = zeros(num_bins,1);
    for bb=1:num_bins
        average_y(bb) = f_avg_y(y(bin_indices{bb}));
    end
end


% ERROR ON Y VALUES
if nargout >=3
    
    % set algorithm
    switch params.error_y
        case 'std_err'
            f_err_y = @(nums)std(nums)/sqrt(length(nums));
        case 'std'
            f_err_y = @(nums)std(nums);
    end

    % fill in variable
    error_y = zeros(num_bins,1);
    for bb=1:num_bins
        error_y(bb) = f_err_y(y(bin_indices{bb}));
    end
end


% NUMBER OF POINTS IN EACH BIN
if nargout >=4

    % fill in variable
    bin_counts = zeros(num_bins,1);
    for bb=1:num_bins
        bin_counts(bb) = length(bin_indices{bb});
    end
end

    
    
    
    