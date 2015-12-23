function [datarun extras] = load_cones(datarun, ind, varargin)
% LOAD_CONES    Find folder with single cone data and load them
% usage: datarun = load_cones(datarun, [ind])
%
% Calls FIND_CONE_DATA to get all matching folders.  If there is only one
% matching folder, loads those data.  
%
% If there are multiple folders and no IND is specified, prints the 
% possible folders and exits. If there are multiple folders and IND is
% specified and numerical, loads the data from the folder index IND.  If
% IND is text, loads the data from the first folder that has a substring
% matching IND (using STRFIND).
% 
% See also FIND_CONE_DATA, STRFIND, IMPORT_SINGLE_CONE_DATA
%
% 2011-07 phli
%
 
% path2data=datarun.names.rrs_prefix;
% tmp=regexp(path2data,'data');
% path2data=path2data(1:tmp(1)-1);

cone_data = find_cone_data(datarun);
if isempty(cone_data)
    warning('No cone data found');
    return;
end

if length(cone_data) > 1 && (nargin < 2 || isempty(ind))
    disp('Multiple matching cone data found:');
    disp_cone_data(cone_data);
    disp('Run again with desired cone data index as 2nd argument.');
    return;
end

if length(cone_data) == 1 && nargin < 2
    datarun.names.cones = cone_data{1};
elseif isnumeric(ind) && ~isempty(ind)
    datarun.names.cones = cone_data{ind};
elseif ischar(ind)
    datarun.names.cones = detect(cone_data, @(cd) (regexp(cd, ind)));
    if isempty(datarun.names.cones)
        disp(['No cone data found matching string "' ind '". Data found:']);
        disp_cone_data(cone_data);
        return;
    end
end
 
[datarun extras] = import_single_cone_data(datarun, datarun.names.cones, varargin{:});
 
% PHL: I don't really see the point of keeping the extras separate.  Move the
% things I care about into datarun
if isfield(extras, 'results') && isfield(extras.results, 'bcf')
    if isfield(extras.results.bcf, 'dll')
        datarun.cones.dll = extras.results.bcf.dll;
    end
end



function disp_cone_data(cone_data)
for i = 1:length(cone_data)
    disp(['  ' num2str(i) ') ' cone_data{i}]);
end