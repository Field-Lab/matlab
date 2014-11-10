function recalc_coord(cone_list)

global new_cone_list
%
% to get coordinates: 
% flip cone image left/right (holistic view, 6.5x?);
% read in small focusing stimulus picture (same magnification), flip
% left/right; get coordinates of the cross, e.g. xcenter=770, ycenter=393
% read in medium fine mapping picture (same magnif); use ginput to measure
% sizes of indiv.checkers (10pxl on the screen). E.g. size is 100 matlab
% pxl - then 1 stim pxl is 10 matlab pxl, magnfactor=10;
% for 320x320 field size (CRT), the whole picture would be 3200x3200 matlab
% pxl
% if a cone has coordinates xcone=600, ycone=305, center (770,393), then stimulus
% location would be: (xcenter-xcone)=170; (ycenter-ycone)=88;
% 170/magnfactor=17; 88/magnfactor=8.8; these are coordinates on the
% stimulus from the center (++ is top left, +- is bottom left, -+ is top right, -- is bottom right)
% if coordinates (0,0) is top left corner, xconeScreen=160-17=143;
% yconeScreen=160-8.8=152;

xcenter=577;
ycenter=303;
magnfactor=10;


xcones=160-round((xcenter-cone_list(:,1))/magnfactor);
ycones=160-round((ycenter-cone_list(:,2))/magnfactor);

new_cone_list=[xcones ycones];
