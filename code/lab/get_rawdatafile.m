function rdf = get_rawdatafile(data_spec)
% For now just uses data_spec as a direct path, but in future would like to
% make this more flexible.  The path can be either to a specific bin file
% or else to a directory of bins where the directory is treated as a single
% "virtual" data file as it is in Vision.

rdf = edu.ucsc.neurobiology.vision.io.RawDataFile(data_spec);