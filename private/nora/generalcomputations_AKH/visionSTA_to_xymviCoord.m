% AKHeitma  2013-12-08
%Only works for same dimsension byt potentially different pixel size
% datarun.vision.sta_fits.mean
% datarun.vision.sd
%
%

function [ center , sd ] = visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, masterdim, slvdim)

x_coord   = round( stafit_centercoord(1)* (slvdim.width  /masterdim.width)  );
y_coord   = slvdim.height - round( stafit_centercoord(2)* (slvdim.height /masterdim.height) );

center.x_coord = x_coord;
center.y_coord = y_coord;

sd.xdir = round( stafit_sd(1)* (slvdim.width   / masterdim.width)  );
sd.ydir = round( stafit_sd(2)* (slvdim.height  / masterdim.height)  );

end


