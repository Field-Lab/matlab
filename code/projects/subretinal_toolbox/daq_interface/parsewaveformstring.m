function [M, param_names] = parsewaveformstring(string)
% [M, param_names] = parsewaveformstring(string)
%
% This function parses a Stimulus.ToString string into two cell arrays that
% can be used to process a dataset, or to load a stimulus from a logfile.
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
param_names = {'Pulse duration',...
               'Pulse times',...
               'Frac. of max. irradiance'};
           
delimiters = strfind(string,';');

% Parsing all the fields in order, starting by 'Start time'. 

% Pulse duration
M{1} = str2num(string(1:delimiters(1)-19));                                 %#ok<ST2NM>

% Pulse time
M{2} = str2num(string(delimiters(1)+15:delimiters(2)-16));                  %#ok<ST2NM>

% Frac of max irradiance
M{3} = str2num(string(delimiters(2)+2:end-24));                             %#ok<ST2NM>

end % parselogfileline