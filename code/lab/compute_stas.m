function stas = compute_stas(spikes, first_ttl, buffered_movie, sta_depth, sta_offset)
% COMPUTE_STAS    Compute many STAs in Vision, convert to Matlab format
%
% usage: stas = compute_stas(spikes, first_ttl, buffered_movie, sta_depth)
%
% inputs: SPIKES            This is a cell array.  Each element holds an 
%                           int32 vector of the spikes for a single STA.  
%                           The spike times are given as integers, so you 
%                           would typically build them as:
%                               spike_times_in_seconds * sampling_rate
%
%         FIRST_TTL         This is the first trigger time, again as an 
%                           integer so multiply the time in seconds by the 
%                           sampling rate
%
%         BUFFERED_MOVIE    This is a Java object of the buffered_movie
%                           from Vision.  You can get this with: 
%                               compute_vision_movie(xml_filename, triggers)
%
%         STA_DEPTH         How many frames long is the STA
%
%         STA_OFFSET        How long after the spike is the last frame.
%                           Defaults to -1, which means the spike is in the
%                           last frame.  0 would give one extra frame after 
%                           the spike at the end.  (May not be accurate for 
%                           interval > 1.)
%
% The Vision code tries to be efficient by only looping through
% BUFFERED_MOVIE once.  Probably could be made more efficient on the Java
% end.
%
% 2010-02 phli
%

if nargin < 5
    sta_offset = -1;
end

import edu.ucsc.neurobiology.vision.matlab.*;
java_stas = Matlab.computeSTAs(spikes, first_ttl, buffered_movie, sta_depth, sta_offset);

stas = cell(length(java_stas), 1);
for i = 1:length(java_stas)
    stas{i} = sta_from_java_sta(java_stas(i));
end