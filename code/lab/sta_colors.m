function color_triplets = sta_colors(sta, sig_stixels, varargin)
% sta_colors     report the color of stixels at significant stixels
%
% usage:  datarun = sta_colors(datarun, sig_stixels, params)
%
% arguments:           sta - 4D array
%              sig_stixels - single indices identifying the significant stixels
%                            
%
% outputs:  color_triplets - N by 3 matrix identifying the RGB triplets
%                            for each of the N signifcant stixels
%
%
% optional arguments            default          explanation 
%           verbose              0               logical for determining optional verbose output  
%
%
% GDF 2008-10-2
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addRequired('sta', @isnumeric);
p.addRequired('sig_stixels', @isnurmic);

p.addParamValue('verbose', false,@islogical)

% parse the inputs
p.parse(sta, sig_stixels, varargin{:});
verbose = p.Results.verbose;

if verbose
    fprintf('\nComputing something important...');
    start_time = clock; % note when it started
end


% ERROR CHECKING    
if length(size(sta)) ~= 4
    error('sta is not a 4D array')
end

if size(sta,3) == 1
    warning('sta is black and white')
end


% BODY OF THE FUNCTION
num_stixels = length(sig_stixels);
temp_color_triplets = zeros(num_indices);

for stixel = 1:num_stixels
    stixel_color = sta(sig_stixel(stixel));
    temp_color_triplets(stixel,:) = stixel_color;
end

color_triplets = temp_color_triplets;  


% display how long it took
if verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end