function cone_data = find_cone_data(datarun)
% FIND_CONE_DATA    Look in SINGLE_CONE_PATH for folder with matching piece and run
% usage: cone_data = find_cone_data(datarun)
%
% Returns a cell array with all matching folder names
%
% 2011-07 phli
%

parsed = parse_rrs_prefix(datarun);

path2data=datarun.names.rrs_prefix;
tmp=regexp(path2data,'data');
path2data=path2data(1:tmp(1)-1);

tmp=regexp(path2data,'streamed');
if ~isempty(tmp); path2data=path2data(1:tmp(1)-1); end

tmp=regexp(path2data,'Streamed');
if ~isempty(tmp); path2data=path2data(1:tmp(1)-1); end

cone_data = find_cone_data(parsed.piece_fullname, parsed.run_name,path2data);