function datarun = stimulus_from_wnmoviefile(datarun)

wnmoviefile = edu.ucsc.neurobiology.vision.io.chunk.WhiteNoiseMovieFile(datarun.names.rrs_movie_path, edu.ucsc.neurobiology.vision.io.chunk.ChunkFile.READ);
stimulus = stimulus_from_wnmoviefile(wnmoviefile);
wnmoviefile.close;

if ~isfield(datarun, 'stimulus') || isempty(datarun.stimulus)
    datarun.stimulus = struct();
end
datarun.stimulus = setstructfields(datarun.stimulus, stimulus);