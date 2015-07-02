function vsta = cell_array_to_vision_sta(carr, r, stixelsz)
%CELL_ARRAY_TO_VISION_STA Converts a cell array STA to the vision format.
%
%  VSTA = CELL_ARRAY_TO_VISION_STA(CARR, R) Converts the sta cell array 
%  CARR to a format that Vision can handle. R is the refresh rate of the
%  STA specified in seconds. Stixel size defaults to 10.
%
%  VSTA = CELL_ARRAY_TO_VISION_STA(..., SZ) Further lets you specify the
%  stixel size.

import edu.ucsc.neurobiology.vision.stimulus.*

if nargin == 2
    stixelsz = 10;
end

vsta = STA(length(carr), size(carr{1},1), size(carr{1}, 2), r, stixelsz, stixelsz);
vstaframe = STAFrame(size(carr{1},1), size(carr{1}, 2), stixelsz, stixelsz);

for k = 1:length(carr)
    cframe = permute(carr{k}, [3, 2, 1]);
    vstaframe = STAFrame(size(carr{1}, 2), size(carr{1}, 1), stixelsz, stixelsz);
    vstaframe.setBuffer(cframe(:) - 0.5);
    vsta.setFrame(k-1, vstaframe);
end

end % cell_array_to_vision_sta

% STAFile(
%         String fileName, int headerCapacity, int width, int height,
%         int staDepth, int staOffset, double stixelWidth, double stixelHeight, double refreshTime)