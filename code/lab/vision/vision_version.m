function [versionint, versionstring] = vision_version()
% VISION_VERSION    Return the version of Vision that is currently on javapath
% usage: [versionint, versionstring] = vision_version()
%
% These version numbers are entered manually into VisionParams.java, so
% they are not guaranteed to be correct.
%
% The versionint will be in the format: 7002000, for example, which 
% translates to string: "7.2.0"
%
% 2010-05 phli
%
import('edu.ucsc.neurobiology.vision.util.*');

versionint = VisionParams.VERSION;

if nargout > 1
    if ismember('versionString', VisionParams.methods)
        versionstring = VisionParams.versionString;
    else
        versionstring = vision_version_to_string(versionint);
    end
end



function versionstring = vision_version_to_string(versionint)
versionstring = [                                 ...
    num2str(floor(versionint / 1000000))          ...
    '.'                                           ...
    num2str(mod(floor(versionint / 1000), 1000))  ...
    '.'                                           ...
    num2str(mod(versionint, 1000))                ...
];
