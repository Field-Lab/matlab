thresh_or_raster=1;
% thresh is 0, raster is 1

previous_streaming_length=50; % in seconds
piece='2014-08-20-0';
data='data002';
number_of_electrodes=512;
raster_interval = 1; % in trigger numbers
raster_length = 2 ; % seconds
number_rasters = 100; % just overestimate - this is for allocation

% Can't actually run these in a row because thresh function ends before you actually
% set the thresh, so make sure thresh is done before you run raster!

if ~thresh_or_raster
    stream_thresh(['/Volumes/stream-bertha/Data/' piece '/' data '.bin'],'thresh.mat',number_of_electrodes)
else
    spikes = stream_raster(['/Volumes/stream-bertha/Data/' piece '/' data '.bin'],'thresh.mat',raster_interval, raster_length, number_of_rasters)
end

% If you want to combine multiple datasets
data2='data003';
spikes_combined = stream_raster(['/Volumes/stream-bertha/Data/' piece '/' data2 '.bin'],'thresh.mat',raster_interval, raster_length, number_of_rasters, spikes)

% If you want an offset
offset = 2; % in number of triggers to skip before doing the rasters
spikes = stream_raster(['/Volumes/stream-bertha/Data/' piece '/' data '.bin'],'thresh.mat',[raster_interval, offset], raster_length, number_of_rasters)
