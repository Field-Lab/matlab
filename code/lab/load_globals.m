function datarun = load_globals(datarun)
% load_globals     load globals file java object in datarun.globals.  overwrites any pre-exisiting java object.
%
% usage:  datarun = load_globals(datarun)
%
% arguments:     datarun - datarun struct
%
% outputs:     datarun - datarun struct with datarun.globals params
%
% 2010-01  gauthier
%

if ~isfield(datarun.names, 'rrs_globals_path') || isempty(datarun.names.rrs_globals_path),
    return
end

% get file name
gf_path = datarun.names.rrs_globals_path;

% if it exists...
if 2 == exist(gf_path,'file')
    % load it
    READ = edu.ucsc.neurobiology.vision.io.chunk.ChunkFile.READ;
    datarun.globals = edu.ucsc.neurobiology.vision.io.chunk.GlobalsFile(gf_path, READ);

    % The current (2012-01) Java implementation of the globals file reads
    % the whole buffer into memory and then closes the file handle already.  
    % Therefore, datarun.globals.close() essentially does nothing, but it 
    % is good practice.
    datarun.globals.close();
end