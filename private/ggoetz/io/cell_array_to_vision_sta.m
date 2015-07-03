function vsta = cell_array_to_vision_sta(ta, e_ta, r, stixelsz)
%CELL_ARRAY_TO_VISION_STA Converts a cell array STA to the vision format.
%
%  VSTA = CELL_ARRAY_TO_VISION_STA(STA, E_STA R) Converts the sta cell 
%  array STA with associated error E_TA to a format that Vision can 
%  handle. R is the refresh rate of the STA specified in seconds. 
%  STA Stixel size defaults to 10.
%  E_STA can be empty.
%
%  VSTA = CELL_ARRAY_TO_VISION_STA(..., SZ) Further lets you specify the
%  stixel size.

import edu.ucsc.neurobiology.vision.stimulus.*

if nargin == 3
    stixelsz = 10;
end
if isempty(e_ta)
    specifyerror = false;
else
    specifyerror = true;
end

vsta = STA(length(ta), size(ta{1},1), size(ta{1}, 2), r, stixelsz, stixelsz);
vstaframe = STAFrame(size(ta{1},1), size(ta{1}, 2), stixelsz, stixelsz);

for k = 1:length(ta)
    cframe = permute(ta{k}, [3, 2, 1]);
    vstaframe = STAFrame(size(ta{1}, 2), size(ta{1}, 1), stixelsz, stixelsz);
    vstaframe.setBuffer(cframe(:));
    if specifyerror
        cerror = permute(e_ta{k}, [3, 2, 1]);
        vstaframe.setErrorBuffer(cerror(:));
    end
    vsta.setFrame(k-1, vstaframe);
end

end % cell_array_to_vision_sta

% STAFile(
%         String fileName, int headerCapacity, int width, int height,
%         int staDepth, int staOffset, double stixelWidth, double stixelHeight, double refreshTime)