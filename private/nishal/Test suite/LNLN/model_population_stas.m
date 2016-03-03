function model=model_population_stas(varargin)

%% Model cell

p = inputParser;
% specify list of optional parameters
ang = pi/3;
p.addParamValue('coneLatticeOrientation', ang);

gridsz=140;
p.addParamValue('gridsz',gridsz);


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

%% Model parameters

% Data says at eccentricity of 5 mm, cone density is between 5000-10000
% cones/mm2 . In our case, each stixel is 5.5 micrometers size, so for an 8
% stixel resolution movie we see generally, each stixel has roughly half-1 of
% cone .. 


% three times the resolution of general moveis used by us .

gridSzX = params.gridsz; % on pixel in 8 stixel resolution means 24 grid pixels
gridSzY = params.gridsz;
stix = 1 ;

% Cone parameters
% on model : 256x256 grid, conegauss sd = 2, cone spacing =12;
nCones_max = 64*64*64;
coneSpacing = 4; % in grid units 
coneLocSd = 0.25; 
coneLatticeOrientation = params.coneLatticeOrientation; % in radians
coneWindowSz = 4; % actually, its half of cone window.
coneGaussSd = 2; % the cone gaussian input profile
% parameters for cone weights?




%% Generate cone mosaic
d = coneSpacing;

xdim = sqrt(nCones_max);
ydim = sqrt(nCones_max);

% Generate lattice vertices
conesX = [];
conesY = [];

% Rotate lattice
theta = coneLatticeOrientation;
R = [cos(theta),-sin(theta);
    sin(theta),cos(theta)];

for ix = -xdim:xdim
    for iy=-ydim:ydim
        conx = (ix-1)*sqrt(3)*d/2+10;
        cony = (d/2)*rem(ix,2) + (iy-1)*d;
        % rotate
        rot = R*[conx;cony];
        rot = rot+randn(2,1)*coneLocSd;
        
        % jitter
        conx=rot(1);cony=rot(2);
        
        if(conx <=gridSzX & cony <=gridSzY & conx>=1 & cony>=1)
        conesX = [conesX;conx];
        conesY = [conesY;cony];
        end
        
    end
end

nCones = size(conesX,1);
 
% conesX = conesX-min(conesX)+10;
% conesY = conesY-min(conesY)+10;

figure;
plot(conesX,conesY,'.');
title('Cone Locations');
axis equal
hold on
n=100;
angle = 0:2*pi/n:2*pi;            % vector of angles at which points are drawn
R = 1*coneGaussSd;                         % Unit radius
x = R*cos(angle);  y = R*sin(angle);   % Coordinates of the circle
for icone=1:nCones
    plot(x+conesX(icone),y+conesY(icone),'r');                      % Plot the circle 
    axis equal;
    grid on;
end

% Make cone input weights
cone_data=cell(nCones,1);
totalConeMap = sparse(gridSzX,gridSzY);

stas = zeros(nCones,numel(totalConeMap));

for icone = 1:nCones
coneMap = zeros(gridSzX,gridSzY);
coneCenterX = round(conesX(icone));
coneCenterY = round(conesY(icone));

coneMap(coneCenterX,coneCenterY) = 1;
coneMap = imgaussfilt(coneMap , coneGaussSd^2);
coneMap=coneMap/max(coneMap(:));

coneMap=coneMap(1:gridSzX,1:gridSzY);

coneWindow = 0*coneMap;
coneWindow(max(coneCenterX-coneWindowSz,1):min(coneCenterX+coneWindowSz,gridSzX),max(coneCenterY-coneWindowSz,1):min(coneCenterY+coneWindowSz,gridSzY))=1;
cone_data{icone}.coneWindow = logical(coneWindow);
cone_data{icone}.coneMap = coneMap;
stas(icone,:) = coneMap(:);
totalConeMap = totalConeMap + coneMap;
end

figure;
imagesc(totalConeMap');
colormap gray
title('Cone map');
% n=100;
% angle = 0:2*pi/n:2*pi;            % vector of angles at which points are drawn
% R = 3*coneGaussSd;                         % Unit radiu
% x = R*cos(angle);  y = R*sin(angle);   % Coordinates of the circle
% for icone=1:nCones
%     hold on;
%     plot(x+conesX(icone),y+conesY(icone),'r');                      % Plot the circle 
%     axis equal;
%     grid on;
% end


% Overlay cone map with stixels
lowRes = gridSzX/stix ;

if(rem(lowRes,2)==0)
lowRes=lowRes+1;
end

% Generate low resolution checker board
checkerboard = rem(reshape(1:lowRes*lowRes,[lowRes,lowRes]),2);
checkerboard=checkerboard(1:end-1,1:end-1);
stimuli = repelem(checkerboard,stix,stix)-0.5;
stimConeMap = stimuli*0.5+totalConeMap;

figure;
imagesc(stimConeMap);
colormap gray
title(sprintf('%d stixel stimulus overlayed with cone map',stix));



%% Temporal filter
Filtlen = 30;
scale_one=1;
scale_two=0.25;
tau_one=4;
tau_two=10;
n_filters=6;
t=[0:Filtlen-1];
ttf = scale_one*((t/tau_one).^n_filters).*exp(-n_filters*(t/tau_one -1)) - scale_two*((t/tau_two).^n_filters).*exp(-n_filters*(t/tau_two -1));
ttf=ttf-mean(ttf);

figure;
plot(ttf)
title('Temporal Filter');

%% Make model strucutre

model.gridSzX = gridSzX;
model.gridSzY = gridSzY;

model.nCones=nCones;
model.cone_data=cone_data;
model.totalConeMap = totalConeMap;
model.conesX = conesX;
model.conesY = conesY;

model.coneSpacing = coneSpacing; % in grid units
model.coneLocSd = coneLocSd; 
model.coneLatticeOrientation = coneLatticeOrientation; % in radians
model.coneWindowSz=coneWindowSz;
model.coneGaussSd = coneGaussSd; % the cone gaussian input profile
% parameters for cone weights?
model.stas = stas;


% temporal filter
model.ttf=ttf';

% scale
model.scale_compared_to_usual_stix = 1;


end
