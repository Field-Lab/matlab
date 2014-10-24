function foam_export_movie_obvius(file, stim, params);
% MKLUMCONSTIM_EXPORT  Generate movie file for experiments
%
% usage: state = mkLumConStim_export(file, stim_params, frames);
%        state = mkLumConStim_export(..., state);
%
% arguments:        file - path of new movie file
%            stim_params - structure of static stimulus paramters
%                 frames - number of stimulus frames
%                  state - random state of algorithm
%
% Export mkLumConStim stimulus (i.e. bubbles) into EJ movie
% file format. Returns state of algoirthm to let you generate new
% stimulus where this left off.
%
% shlens 2009-02-24
%

% frames
frames = params.dims(3);

% generate header
version = 1.1;%foam('version');
header = generate_foam_header(params, version);
fprintf('(foam_export_movie_obvius) writing frames\n');

% start file
fid = fopen(file, 'w');
fprintf(fid, header);
T=text_waitbar;
% do it
for i =1:frames
   T=text_waitbar(T,i/frames - 0.01);
  % reshape matrix so monitor rows are successfully written
  % Note: default vectorizing in Matlab is column-wise
  stim_frame = reshape(stim(:,:,i)',[1 numel(stim(:,:,i))]); 

  % write it out in EJ lab format as RGB triples (grey scale)
  fwrite(fid, repmat(256*stim_frame, [3 1]), 'ubit8');

end
fprintf('\n');
% close file
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hdr = generate_foam_header(params, version)

hdr = sprintf('width\t%d\r\n', params.dims(2));
hdr = [hdr, sprintf('height\t%d\r\n', params.dims(1))];
hdr = [hdr, sprintf('frames-generated\t%d\r\n', params.dims(3))];
if params.contrast_sd == 1
  hdr = [hdr, sprintf('algorithm\tcrest\r\n')];
else
  hdr = [hdr, sprintf('algorithm\tfoam\r\n')];
end
hdr = [hdr, sprintf('version %.2f, %s\r\n', version, date)];
hdr = [hdr, sprintf('randstate\t%d\r\n', params.seed)];
hdr = [hdr, sprintf(['dur_sec\t ' mat2str(params.dur_sec,4) '\r\n'])];
hdr = [hdr, sprintf(['gauss_spatial_sd\t ' mat2str(params.gauss_spatial_sd,4) '\r\n'])];
hdr = [hdr, sprintf(['tau_temp_interval\t ' mat2str(params.tau_temp_interval,4) '\r\n'])];
if params.contrast_sd == 1
  hdr = [hdr, sprintf(['thresholding\t ' params.thresholding '\r\n'])];
else
  hdr = [hdr, sprintf(['contrast_sd\t ' mat2str(params.contrast_sd,4) '\r\n'])];
end
hdr = [hdr, sprintf(['\r\n'])];
