function model=model_LNLN_fewSU()
%% Model cell

%% Model parameters

% Data says at eccentricity of 5 mm, cone density is between 5000-10000
% cones/mm2 . In our case, each stixel is 5.5 micrometers size, so for an 8
% stixel resolution movie we see generally, each stixel has roughly half-1 of
% cone .. 


% three times the resolution of general moveis used by us .

gridSzX = 480; % 8 stixel resolution means 24 grid pixels
gridSzY = 480;
stix = 8 * 3;

% Cone parameters

nCones = 25;
coneSpacing = 30; % in grid units
coneLocSd = 1.5; 
coneLatticeOrientation = pi/3; % in radians
coneWindowSz = 4; % actually, its half of cone window.
coneGaussSd = 2; % the cone gaussian input profile
% parameters for cone weights?


% Sub-unit parameters
nSU = 4; % ~10 in center and 50 in surround
avgConeperSU = 5; % should be 5-10 

% Cone to SU weights 
cone_to_SU_Wt = 1*ones(nSU,nCones)/(45);

% sample from bipolar to ganglion weights from uniform distribution.
% SU_gang_weights = 0.5+(rand(nSU,1)/2);
sdSUweights = 2; % ideally 0.5

% Sub-unit NL
fs = @(x) exp(x);

% ganglion cell NL
g = @(x) x - min(x);

%% Generate cone mosaic
d = coneSpacing;

xdim = sqrt(nCones);
ydim = sqrt(nCones);

% Generate lattice vertices
conesX = [];
conesY = [];
for ix = 1:xdim
    for iy=1:ydim
    conesX = [conesX;(ix-1)*sqrt(3)*d/2+10];
    conesY = [conesY;(d/2)*rem(ix,2) + (iy-1)*d];
    end
end

% Rotate lattice
theta = coneLatticeOrientation;
R = [cos(theta),-sin(theta);
    sin(theta),cos(theta)];

A = R*[conesX,conesY]';
conesX = A(1,:)';
conesY = A(2,:)';

% Add jitter 
conesX = conesX + randn(nCones,1)*coneLocSd;
conesY = conesY + randn(nCones,1)*coneLocSd;

conesX = conesX-min(conesX)+10;
conesY = conesY-min(conesY)+10;

figure;
plot(conesX,conesY,'*');
title('Cone Locations');
axis equal

% Make cone input weights
cone_data=cell(nCones,1);
totalConeMap = zeros(gridSzX,gridSzY);

for icone = 1:nCones
coneMap = zeros(gridSzX,gridSzY);
coneCenterX = round(conesX(icone));
coneCenterY = round(conesY(icone));
coneMap(coneCenterX,coneCenterY) = 1;
coneMap = imgaussfilt(coneMap , coneGaussSd^2);
coneMap=coneMap/max(coneMap(:));

coneWindow = 0*coneMap;
coneWindow(max(coneCenterX-coneWindowSz,1):min(coneCenterX+coneWindowSz,gridSzX),max(coneCenterY-coneWindowSz,1):min(coneCenterY+coneWindowSz,gridSzY))=1;
cone_data{icone}.coneWindow = logical(coneWindow);
cone_data{icone}.coneMap = coneMap;
totalConeMap = totalConeMap + coneMap;
end

figure;
imagesc(totalConeMap);
colormap gray
title('Cone map');

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
imagesc(stimConeMap(1:floor(max(conesX))+10,1:floor(max(conesY))+10));
colormap gray
title(sprintf('%d stixel stimulus overlayed with cone map',stix));


%% Sub-units 

% Compute matrix for going from cones to sub-units 

[coneidx,C] = kmeans([conesX,conesY],nSU);
cols = distinguishable_colors(nSU+1);

figure;
for isu=1:nSU
plot(conesX(coneidx==isu),conesY(coneidx==isu),'*','Color',cols(isu,:));
hold on;
end

connectionMat = zeros(nCones,nSU);
for isu = 1:nSU
connectionMat(:,isu) = (coneidx==isu);
end

%% Cone to SU weights
connectionMat = connectionMat';
connectionMat = cone_to_SU_Wt.*connectionMat;

cols = distinguishable_colors(nSU,[0.2,0.2,0.2]);
cols = cols(randperm(size(cols,1)),:);
% Make totalConeMap3D 
totalConeMap3D = zeros(gridSzX,gridSzY,3);

listC = 1:nCones;
for isu = 1:nSU
    listC(coneidx==isu)
    for icone = listC(coneidx==isu)
        for idim=1:3
    totalConeMap3D(:,:,idim) = totalConeMap3D(:,:,idim) + double(cone_data{icone}.coneMap>0.2)*cols(isu,idim);
        end
    end
end
%totalConeMap3D = totalConeMap3D*10;

figure;
imagesc(totalConeMap3D);
%% SU to ganglion cell weights
% taken care of in SU_gang_weightsMat

SU_gang_weights = 0.05* 1./sqrt(sum(connectionMat,2));
SU_gang_weights = SU_gang_weights + rand(length(SU_gang_weights),1)*sdSUweights;

%% Non linearities ?

% taken care of by fs and g

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
model.totalConeMap3D = totalConeMap3D;
model.conesX = conesX;
model.conesY = conesY;

model.coneSpacing = coneSpacing; % in grid units
model.coneLocSd = coneLocSd; 
model.coneLatticeOrientation = coneLatticeOrientation; % in radians
model.coneWindowSz=coneWindowSz;
model.coneGaussSd = coneGaussSd; % the cone gaussian input profile
% parameters for cone weights?


% Sub-unit parameters
model.nSU = nSU; 
model.avgConeperSU = avgConeperSU;

% sample from bipolar to ganglion weights from uniform distribution.
model.SU_gang_weights = SU_gang_weights;

% Sub-unit NL
model.fs = fs;

% ganglion cell NL
model.g = g;

model.cone_to_SU_connection = connectionMat; 
model.cone_su_idx=coneidx;


% temporal filter
model.ttf=ttf';

% scale
model.scale_compared_to_usual_stix = 3;


end
