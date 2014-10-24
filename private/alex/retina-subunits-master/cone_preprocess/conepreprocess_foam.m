function datarun = conepreprocess_foam(datarun, varargin)
% CONEPREPROCESS_FOAM    Setup datarun.stimulus for CONE_PREPROCESS on a foam/crest movie
% usage: datarun = conepreprocess_foam(datarun, opts)
%
% Sets up DATARUN.STIMULUS by loading the java_movie and then setting a 
% few additional fields.  OPTS are simply passed on to CONE_PREPROCESS.
%
% See also: CONE_PREPROCESS, CALC_CONE_INPUTS, GET_BINNED_SPIKES
%
% 2012-09 freeman
%

% manually set parameters for the data run
datarun.stimulus.duration = 1800;
datarun.stimulus.refresh_time = 5/60.35; % SET
datarun.stimulus.field_height = 300; % SET
datarun.stimulus.field_width = 300; % SET

% properties of the movie
datarun.stimulus.movie.field_height = datarun.stimulus.field_height;
datarun.stimulus.movie.field_width = datarun.stimulus.field_width;
datarun.stimulus.movie.path = ...
    '/share/users-archive/freeman/crest_b/crest_s1_t2_b_d1800.rawMovie'; % SET
datarun.stimulus.stimobject = datarun.stimulus.movie;
datarun.stimulus.get_movie_frame = @get_movie_frame_foam;
datarun.stimulus.wc_preprocess = @wc_preprocess_javamovie_rawmovie;

% do the cone preprocessing
datarun = cone_preprocess(datarun, varargin{:});

function frame = get_movie_frame_foam(movie, i)

fid = fopen(movie.path);
% parse through the header
fseek(fid,0,'bof');
tline = fgetl(fid);
while ~strcmp(tline,'')
  tline = fgetl(fid);
end
% seek to the desired frame
fseek(fid,(i-1)*(movie.field_width*movie.field_height*3),'cof');
% grabe the frame
stim_frame = fread(fid,[movie.field_width*movie.field_height*3]);
stim_frame = stim_frame / 256;

% do reshaping and traposing to invert the manipulation done
% when writing the file to OBVIUS format... but we don't need to do this
% because we're handling it in the manipulations to the cone matrix 
% leaving it here commented out for now 
% stim_frame = reshape(stim_frame,[3 movie.field_width movie.field_height]);
% stim_frame = squeeze(stim_frame(1,:,:))';
% stim_frame = repmat(stim_frame,[1 1 3]);

% vectorize
frame = stim_frame(:);
% close the file
fclose(fid);