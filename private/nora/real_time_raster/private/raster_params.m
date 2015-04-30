% TRIGGERSv
% In a continuous run, triggers=seconds*120/100

% raster_interval: the triggers at which to start
% raster_length: in seconds
% number_of_rasters: just to preallocate the variable, so overestimate

% In triggers

% from time
dist = 135;
offset = 363;
raster_start = [1:dist:(dist*4) offset+(1:dist:(dist*4))];

% In seconds
raster_length = 6; 
