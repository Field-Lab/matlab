function [xCoords, yCoords] = getElectrodeCoords519()

% coordinate of each electrode on the 519-electrode array, in microns with
% the origin at the center of the array

 xCoords = [330;210;240;360;300;270;330;150;210;240;360;300;120;120;270;330;180;210;240;360;300;90;150;270;330;180;210;240;360;300;90;150;270;330;180;210;240;360;300;120;150;270;330;180;210;240;360;300;60;30;270;330;120;180;240;360;300;90;150;270;330;210;120;240;360;300;180;210;330;270;150;90;300;360;240;180;120;330;270;60;60;300;360;240;210;180;330;270;150;120;300;360;240;210;180;330;270;150;90;300;360;240;210;180;330;270;150;90;300;360;240;210;180;330;270;120;120;300;360;240;210;150;330;270;300;360;240;210;330;300;270;270;240;240;210;210;210;180;180;180;180;180;150;150;150;150;150;150;120;120;120;120;120;120;120;60;90;90;90;90;90;90;90;90;90;60;60;60;60;60;60;60;60;60;60;30;30;30;30;30;30;30;30;30;30;30;30;30;1.83697009358929e-15;5.51091070428436e-15;9.18485089146295e-15;1.28587910786415e-14;1.65327312658201e-14;2.02066723000317e-14;2.38806133342432e-14;2.20436428171374e-14;1.83697017829259e-14;1.46957607487144e-14;1.10218214085687e-14;7.34788037435718e-15;3.67394018717859e-15;-30;-30;-30;-30;-30;-30;-30;-30;-30;-30;-30;-30;-30;-60;-60;-60;-60;-60;-60;-60;-60;-60;-60;-90;-90;-90;-90;-90;-90;-90;-90;-90;-60;-120;-120;-120;-120;-120;-120;-120;-150;-150;-150;-150;-150;-150;-180;-180;-180;-180;-180;-210;-210;-210;-240;-240;-270;-300;-270;-330;-210;-240;-360;-300;-270;-330;-150;-210;-240;-360;-300;-120;-120;-270;-330;-180;-210;-240;-360;-300;-90;-150;-270;-330;-180;-210;-240;-360;-300;-90;-150;-270;-330;-180;-210;-240;-360;-300;-120;-150;-270;-330;-180;-210;-240;-360;-300;-60;-60;-270;-330;-120;-180;-240;-360;-300;-90;-150;-270;-330;-210;-180;-300;-360;-240;-120;-210;-330;-270;-150;-90;-300;-360;-240;-180;-120;-330;-270;-30;-60;-300;-360;-240;-210;-180;-330;-270;-150;-120;-300;-360;-240;-210;-180;-330;-270;-150;-90;-300;-360;-240;-210;-180;-330;-270;-150;-90;-300;-360;-240;-210;-180;-330;-270;-120;-120;-300;-360;-240;-210;-150;-330;-270;-300;-360;-240;-210;-330;-300;-270;-270;-240;-240;-210;-210;-210;-180;-180;-180;-180;-180;-150;-150;-150;-150;-150;-150;-120;-120;-120;-120;-120;-120;-120;-60;-90;-90;-90;-90;-90;-90;-90;-90;-90;-60;-60;-60;-60;-60;-60;-60;-60;-60;-60;-30;-30;-30;-30;-30;-30;-30;-30;-30;-30;-30;-30;0;-3.67394018717859e-15;-7.34788037435718e-15;-1.10218214085687e-14;-1.46957607487144e-14;-1.83697017829259e-14;-2.20436428171374e-14;-2.38806133342432e-14;-2.02066723000317e-14;-1.65327312658201e-14;-1.28587910786415e-14;-9.18485089146295e-15;-5.51091070428436e-15;-1.83697009358929e-15;30;30;30;30;30;30;30;30;30;30;30;30;60;60;60;60;60;60;60;60;60;60;90;90;90;90;90;90;90;90;90;60;120;120;120;120;120;120;120;150;150;150;150;150;150;180;180;180;180;180;210;210;210;240;240;270;270;300];
 yCoords = [-225;-195;-210;-210;-210;-195;-195;-135;-165;-180;-180;-180;-120;-90;-165;-165;-150;-135;-150;-150;-150;-75;-105;-135;-135;-120;-105;-120;-120;-120;-45;-75;-105;-105;-90;-75;-90;-90;-90;-60;-45;-75;-75;-60;-45;-60;-60;-60;-30;-15;-45;-45;-30;-30;-30;-30;-30;-15;-15;-15;-15;-15;-7.34788037435718e-15;-1.46957607487144e-14;-2.20436428171374e-14;-1.83697017829259e-14;-1.10218214085687e-14;15;15;15;15;15;30;30;30;30;30;45;45;-3.67394018717859e-15;30;60;60;60;45;60;75;75;45;60;90;90;90;75;90;105;105;75;45;120;120;120;105;120;135;135;105;75;150;150;150;135;150;165;165;90;120;180;180;180;165;135;195;195;210;210;210;195;225;240;225;255;270;240;285;255;225;240;270;300;210;180;165;195;225;255;285;315;150;180;210;240;270;300;330;60;105;135;165;195;225;255;285;315;345;90;120;150;180;210;240;270;300;330;360;105;135;165;195;225;255;285;315;345;375;75;45;15;30;90;150;210;270;330;390;360;300;240;180;120;60;15;45;75;375;345;315;285;255;225;195;165;135;105;360;330;300;270;240;210;180;150;120;90;345;315;285;255;225;195;165;135;105;60;330;300;270;240;210;180;150;315;285;255;225;195;165;180;210;300;270;240;225;255;285;240;270;255;240;225;225;195;210;210;210;195;195;135;165;180;180;180;120;90;165;165;150;135;150;150;150;75;105;135;135;120;105;120;120;120;45;75;105;105;90;75;90;90;90;60;45;75;75;60;45;60;60;60;30;3.67394018717859e-15;45;45;30;30;30;30;30;15;15;15;15;15;1.10218214085687e-14;1.83697017829259e-14;2.20436428171374e-14;1.46957607487144e-14;7.34788037435718e-15;-15;-15;-15;-15;-15;-30;-30;-30;-30;-30;-45;-45;-15;-30;-60;-60;-60;-45;-60;-75;-75;-45;-60;-90;-90;-90;-75;-90;-105;-105;-75;-45;-120;-120;-120;-105;-120;-135;-135;-105;-75;-150;-150;-150;-135;-150;-165;-165;-90;-120;-180;-180;-180;-165;-135;-195;-195;-210;-210;-210;-195;-225;-240;-255;-225;-270;-240;-285;-255;-225;-240;-270;-300;-210;-180;-165;-195;-225;-255;-285;-315;-150;-180;-210;-240;-270;-300;-330;-60;-105;-135;-165;-195;-225;-255;-285;-315;-345;-90;-120;-150;-180;-210;-240;-270;-300;-330;-360;-105;-135;-165;-195;-225;-255;-285;-315;-345;-375;-75;-45;0;-60;-120;-180;-240;-300;-360;-390;-330;-270;-210;-150;-90;-30;-45;-75;-375;-345;-315;-285;-255;-225;-195;-165;-135;-105;-360;-330;-300;-270;-240;-210;-180;-150;-120;-90;-345;-315;-285;-255;-225;-195;-165;-135;-105;-60;-330;-300;-270;-240;-210;-180;-150;-315;-285;-255;-225;-195;-165;-180;-210;-300;-270;-240;-225;-255;-285;-240;-270;-225;-255;-240];