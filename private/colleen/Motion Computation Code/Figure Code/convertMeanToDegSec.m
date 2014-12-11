function mean = convertMeanToDegSec(mean)
% 5 um/pix
% 1deg/225 um
% 10 pix/1 stx (everything converted to this ratio before analysis)
mean = mean*5/225*10;