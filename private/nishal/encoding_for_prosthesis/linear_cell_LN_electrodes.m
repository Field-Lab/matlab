%% Linear cells, LN stimulation model

%%
addpath(genpath('../Test suite/'));
%%
% make on and off parasol mosaic
spatial_extent=64;
on = model_population_stas('coneLatticeOrientation',0/3,'gridsz',spatial_extent);
off = model_population_stas('coneLatticeOrientation',pi/6,'gridsz',spatial_extent);

% plot cells
n=100;
angle = 0:2*pi/n:2*pi;            % vector of angles at which points are drawn
R = 3*on.coneGaussSd;                         % Unit radius
x = R*cos(angle);  y = R*sin(angle);   % Coordinates of the circle

figure;
%plot(on.conesX,on.conesY,'r*');
hold on;

for icone=1:on.nCones
    plot(x+on.conesX(icone),y+on.conesY(icone),'r');                      % Plot the circle 
    axis equal;
    grid on;
end
hold on;

R = 3*off.coneGaussSd;                         % Unit radius
x = R*cos(angle);  y = R*sin(angle);   % Coordinates of the circle

%plot(off.conesX,off.conesY,'b*');
hold on;
n=100;
for icone=1:off.nCones
    plot(x+off.conesX(icone),y+off.conesY(icone),'b');                      % Plot the circle 
    axis equal;
    grid on;
end

%% electrode map 

elecSpacing=4;
arrSz=round(spatial_extent*1.7/elecSpacing);
elecLatticeOrientation = pi/6;
elecs = getElectrodes_simulation(elecSpacing,arrSz,elecLatticeOrientation)
nElecs = length(elecs.x);


hold on;
plot(elecs.x,elecs.y,'.')

on = map_population_electrodes(on,elecs);
off = map_population_electrodes(off,elecs);

figure;
for icell=1:30
subplot(6,5,icell)
plot_cells_elecs(on,elecs,icell)
set(gca,'xTick',[]);set(gca,'yTick',[]);
end

%% Experimentation

gridSzX =  on.gridSzX;
gridSzY = on.gridSzY;

stas_on = on.stas; 
stas_off = off.stas;
stas = [stas_on;-stas_off];
stas_inv = pinv(stas);
nCells = size(stas,1);


weight_elecs = [on.elecs.weight_elecs;off.elecs.weight_elecs];
nl_ec = @(x) 1./(1+exp(-x));
nl_ec_deri = @(x) (exp(-x)./(1+exp(-x)).^2);
%% Input image
img = imread('/Volumes/Lab/Users/bhaishahster/SIPI database/misc/4.2.04.tiff'); 
img = double(img(:,:,1))/255 - 0.5;
%img = double(img(200:199+gridSzX,200:199+gridSzY));
img = imresize(img,[gridSzX,gridSzY]);
img_flat= img(:);

cell_resp= stas*img_flat;
a = min(cell_resp);
normalized_cell_resp = cell_resp-a;
b = max(normalized_cell_resp);
normalized_cell_resp = normalized_cell_resp/b;
current_old = zeros(nElecs,1);

for iter=1:5
    iter
cvx_begin
variables cell_r_norm(nCells) current(nElecs)
minimize (sum_square(stas_inv*(cell_r_norm*b+a) - img_flat))
subject to 
   nl_ec(weight_elecs*current_old) + nl_ec_deri(weight_elecs*current_old).*(weight_elecs*(current-current_old))== cell_r_norm
   %cell_r_norm<=1
   %cell_r_norm>=0
cvx_end
current_old=current;
end

cell_r_norm = nl_ec(weight_elecs*current);
stim_img = stas_inv*(cell_r_norm*b+a);
stim_img =reshape(stim_img,gridSzX,gridSzY);
figure; 
imagesc(stim_img);axis image;
colormap gray

%% Perfect stimulation

stim_img_perfect = stas_inv*(stas*(img_flat));
stim_img_perfect =reshape(stim_img_perfect,gridSzX,gridSzY);
figure;
imagesc(stim_img_perfect);axis image
colormap gray;
