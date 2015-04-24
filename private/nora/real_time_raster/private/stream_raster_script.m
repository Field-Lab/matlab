thresh_or_raster=1;
% thresh is 0, raster is 1

previous_streaming_length=50; % in seconds
piece='9999-99-99-9';
data='data014';
number_of_electrodes=512;
raster_starts=1:5:200; % in trigger number
raster_length=2; % seconds

% Can't actually run these in a row because thresh function ends before you actually
% set the thresh, so make sure thresh is done before you run raster

if ~thresh_or_raster
    stream_thresh(['/Volumes/stream-bertha/Data/' piece '/' data '.bin'],'thresh.mat',number_of_electrodes)
else
    stream_raster(['/Volumes/stream-bertha/Data/' piece '/' data '.bin'],'thresh.mat',raster_starts,raster_length,previous_streaming_length)
end