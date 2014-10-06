function sta_frame = sta_frame_from_java_sta_frame(java_sta_frame, height, width)
% STA_FRAME_FROM_JAVA_STA_FRAME    Convert from Java STAFrame object into standard SNL-E MatLab STA frame matrix
%
% usage: sta_frame = sta_frame_from_java_sta_frame(java_sta_frame, [height, width])
%
% Jeff figured out the proper conversion from java_sta_frame.getBuffer into
% our format; I abstracted this from sta_from_java_sta and cleaned it up.
%
% 2010-02 phli
%

if nargin < 3
    width = java_sta_frame.getWidth;
end
if nargin < 2
    height = java_sta_frame.getHeight;
end

sta_frame = permute(reshape(java_sta_frame.getBuffer, 3, width, height), [3 2 1]);