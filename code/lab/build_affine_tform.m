function outtform = build_affine_tform(tform)
% BUILD_AFFINE_TFORM    Build an affine transformation from various inputs
%
% usage: outtform = build_affine_tform(tform)
%
% TFORM can be a 2x2 rotation matrix, a 3x3 affine transformation (in which
% case we return it untouched), a scalar rotation in clockwise degrees, or
% a keyword 'xreflect' or 'yreflect'.
%
% 2010-02 phli
%


if ischar(tform) 
    
    if strcmp(tform, 'xreflect')
        outtform = [-1 0 0; 0 1 0; 0 0 1];
        
    elseif strcmp(tform, 'yreflect')
        outtform = [1 0 0; 0 -1 0; 0 0 1];
    end
    
elseif isnumeric(tform)
    
    if numel(tform) == 1
        % Assume rotation in degrees
        rads = tform / 180 * pi;
        outtform = [cos(rads) -sin(rads) 0; sin(rads) cos(rads) 0; 0 0 1];
        
    elseif all(size(tform) == [2 2])
        % Assume rotation matrix
        outtform = vertcat([tform [0; 0]], [0 0 1]);
        
    elseif all(size(tform) == [3 3])
        outtform = tform;
    end
    
end