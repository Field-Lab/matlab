function [M, param_names] = parsestimulusstringlabview(string)
% [M, param_names] = parsestimulusstringlabview(string)
%
% This function parses the first part of a Stimulus.ToString string into 
% two cell arrays that can be used to process a dataset, or to load 
% a stimulus from a logfile.
% It is spaghetti code, but dealing with string in Matlab is messy... 
% It could be worth switching to xml stimuli, instead of using strings to
% summarize things.
%
% Parameters:
%   - string: the string to parse
%
% Returns:
%   - M: cell array with the relevant data parsed 
%   - param_names: name of each of the cells of M

% Version: 0.1 - 2012/09/18
% Author: Georges Goetz, Stanford University
% 

M = {};
param_names = {'Number of trials', 'Frequency', 'Wavelength', 'Start time', 'Duration'};

delimiters = strfind(string, ';');

% Number of trials
M{1} = str2double(string(1:(delimiters(1)-8)));

% Frequency
M{2} = str2double(string(delimiters(1)+2:delimiters(2)-4));

% Wavelength
M{3} = string(delimiters(2)+2:delimiters(3)-1);

% Start time
date_str = string(delimiters(3)+13:delimiters(4)-1);
if ~isempty(date_str)
    M{4} = datenum(date_str,'dd-mmm-yyyy HH:MM:SS:.FFF');
else
    M{4} = 0;
end

% Duration
M{5} = str2double(string(delimiters(4)+9:end-2));

end % parsestimulusstring