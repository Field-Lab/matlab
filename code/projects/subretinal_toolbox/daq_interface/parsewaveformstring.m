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
               'Irradiance'};
           
delimiters = strfind(string,';');
waveformregexp = '\[(.*?)\]';

% Parsing all the fields in order, starting by 'Start time'. 

% Pulse duration
[startind, endind] = regexp(string(1:delimiters(1)), waveformregexp);
M{1} = str2num(string(startind:endind));                                   %#ok<ST2NM>

% Pulse time
[startind, endind] = regexp(string(delimiters(1):delimiters(2)), waveformregexp);
M{2} = str2num(string((startind:endind) + delimiters(1) - 1));             %#ok<ST2NM>

% Frac of max irradiance
[startind, endind] = regexp(string(delimiters(2):end), waveformregexp);
M{3} = str2num(string((startind:endind) + delimiters(2) - 1));             %#ok<ST2NM>

end % parselogfileline