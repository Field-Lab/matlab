function vislocal()
% VISSTABLE     Swap out stable or volatile Vision path for local path
%
% You must have a vision_path_local.m defined somewhere on your path.
% 
% Example vision_path_local.m:
%
%   function paths = vision_path_local()
%   paths = {
%       '/Users/peterli/priv/vision7/Vision.jar'
%       '/Users/peterli/priv/cellfinder/Cell-Finder.jar'
%   };
%

javaswappath(vision_path_local, [vision_path_stable vision_path_volatile]);