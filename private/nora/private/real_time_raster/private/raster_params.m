% TRIGGERS
% In a continuous run, triggers=seconds*120/100

% from time
% dist = 135;
% offset = 363;
% raster_start = [1:dist:(dist*4) offset+(1:dist:(dist*4))];

% 
triggers_per_WNtrial = 36 + 12 + 2;
triggers_per_NSEMtrial = 72 + 24 + 2; 
triggers_per_grey = 12; 

% in triggers
raster_start = 513:triggers_per_NSEMtrial:(515+triggers_per_NSEMtrial*10);
Wraster_start = 1:triggers_per_WNtrial:(triggers_per_WNtrial*10+1);

% In seconds
raster_length = 10; 
