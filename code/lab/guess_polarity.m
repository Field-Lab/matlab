function polarity = guess_polarity(timecourse, varargin)
% GUESS_POLARITY     guess whether a cell is ON or OFF
%
% usage:  polarity = guess_polarity(timecourse, <params>)
%
% arguments:   timecourse - FxC matrix of the timecourse, F = frames, C = colors
%                  params - struct or list of optional parameters (see below)
%
% outputs:       polarity - guessed polarity
%                               +1  ON cell
%                               -1  OFF cell
%                               0   timecourse is all zeros or empty
%
%
% optional fields in params, their default values, and what they specify:
%
% ratio         2       threshold ratio for declaring that one peak is "much" larger than the other.
%
%
% This function guesses whether the timecourse comes from an ON or OFF cell by
% looking at the sizes and times of largest + and largest - peaks.
% If one peak is "much" larger (see above), take that as the polarity.
% If the + and - peaks are comparable in size, assume whichever comes later is 
% the dominant one.
%
% 2008-03  gauthier 
% 2010-03  gauthier, added special conditions for returning early
%
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('ratio', 2);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% return 0 if problems
if isempty(timecourse) || all(reshape(timecourse,[],1)==0) || any(isnan(reshape(timecourse,[],1)))
    polarity = 0;
    return
end



% get timecourse for polarity identification by summing individual timecourses
tcp = sum(timecourse,2);

% get peaks
big_on = max(tcp);
big_off = min(tcp);

% if timecourse doesn't go both positive and negative, use whichever side it's on
if big_on <= 0
    polarity = -1;
    return
else
    if big_off >= 0
        polarity = 1;
        return
    end
end

% if one is larger by a factor of 2
if max(big_on,-big_off)/min(big_on,-big_off) > params.ratio
    % assume it is the polarity of the cell
    polarity = sign(big_on + big_off);
else
    % otherwise use whichever comes later

    % find when the timecourse attained its min and max
    on_times = find(tcp == big_on);
    off_times = find(tcp == big_off);
    
    % determine whether the min or max came later
    % if the timecourse attains an extremum multiple times, use the latest occurence
    polarity = sign(on_times(end) - off_times(end));
end

