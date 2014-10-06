function the_path = server_data_path
% return path to the server where data files are kept
%
% This is for the old SNL-E /snle/data stuff, which is mostly confocal
% images now kept under /Analysis.  We should probably make a separate one
% for server_rawdata_path if that's ever needed.

the_path = '/Analysis/';