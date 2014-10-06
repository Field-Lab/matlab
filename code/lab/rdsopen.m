function rds = rdsopen(filepath, header, varargin)
% RDSOPEN   Open Vision RawDataSaver object
% usage: rds = rdsopen(filepath, header, varargin)
%
% Make sure to call RDSCLOSE(RDS) when you're done! 
%
% 2011-08, phli
%

opts = inputParser();
opts.addParamValue('samples_to_buffer', 100);
opts.addParamValue('num_buffers', 250);
opts.addParamValue('commonpath', ''); % Only useful if FILEPATH refers to multiple files
opts.parse(varargin{:});
opts = opts.Results;

seconds_to_stream = 0; % Write all data; don't stop after given number of seconds.  Don't think this is useful here.

% Open RawDataSaver
rds = edu.ucsc.neurobiology.vision.io.RawDataSaver(filepath, opts.commonpath, header, opts.samples_to_buffer, opts.num_buffers, seconds_to_stream);
