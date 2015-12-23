function cone_data = find_cone_data(piece, run, path2data)
% FIND_CONE_DATA    Look in SINGLE_CONE_PATH for folder whose name includes PIECE _ RUN
% usage: cone_data = find_cone_data(piece, run)
%
% Returns a cell array with all matching folder names
%
% 2011-07 phli
%

% files = dir(single_cone_path());
files = dir(path2data); 

 
cone_data = {};
for i = 1:length(files)
    file = files(i);
    if (~file.isdir), continue; end
    
    if ~isempty(strfind(file.name, [piece '_' run])) || ...
       ~isempty(strfind(file.name, [piece '_streamed_' run])) ...
       ~isempty(strfind(file.name, [piece '_Streamed_' run]))
            cone_data{end+1} = file.name;
    end 
end
