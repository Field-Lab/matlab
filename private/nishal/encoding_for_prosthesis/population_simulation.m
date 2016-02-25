%%
addpath(genpath('../Test suite/'));
%%
% make on and off parasol mosaic

on = model_population_stas('coneLatticeOrientation',0/3,'gridsz',64);
off = model_population_stas('coneLatticeOrientation',pi/6,'gridsz',64);

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

%%
% stas
gridSzX =  on.gridSzX;
gridSzY = on.gridSzY;

stas_on = on.stas; 
stas_off = off.stas;
stas = [stas_on;-stas_off];
nCells = size(stas,1);

% make co-stimulation matrix
% locations
cellsX = [on.conesX;off.conesX];
cellsY = [on.conesY;off.conesY];
cones_loc = [cellsX,cellsY];

co_stim_mat = [];
prob_drop_rate=10;
for icone =1:nCells
    for jcone = 1:nCells
        
        if(icone~=jcone)
    d = norm(cones_loc(icone,:) - cones_loc(jcone,:));
    p = exp(-d^2/prob_drop_rate);

    xx = sparse(nCells,1);
    if(rand()<p)
    xx(icone)=1;
    xx(jcone)=1;
    co_stim_mat = [co_stim_mat,xx]; % add individual cell stimulation ? 
    end
        end
    end
end


if(isempty(co_stim_mat))
co_stim_mat = eye(nCells,nCells);
else
for icone =1:nCells
if(sum(co_stim_mat(icone,:),2)==0)
xx = zeros(nCells,1);
xx(icone) = 1;
co_stim_mat = [co_stim_mat,xx];
end
end
end


figure;
plot(cellsX,cellsY,'k*');hold on;
ii=1:nCells;
for i=1:size(co_stim_mat,2)
xx = logical(co_stim_mat(:,i));
cells = ii(xx);
plot(cones_loc(cells,1),cones_loc(cells,2),'-*');
hold on;
end

stas_inv = pinv(stas);
co_stim_mat_inv = pinv(full(co_stim_mat));

co_stim_mat_perfect = eye(nCells);
co_stim_mat_inv_perfect = pinv(co_stim_mat_perfect);
%% input stimulus ? 
% stim = load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/testmovie_schemeA_8pix_Identity_8pix.mat');
% itime=500;
% img = stim.testmovie.matrix(:,:,itime);

img = imread('/Volumes/Lab/Users/bhaishahster/SIPI database/misc/4.2.04.tiff'); 
img = double(img(:,:,1))/255 - 0.5;
%img = double(img(200:199+gridSzX,200:199+gridSzY));
img = imresize(img,[gridSzX,gridSzY]);

stim_img = stas_inv*(co_stim_mat*(co_stim_mat_inv*stas*(img(:))));
stim_img =reshape(stim_img,gridSzX,gridSzY);

stim_img_perfect = stas_inv*(co_stim_mat_perfect*(co_stim_mat_inv_perfect*stas*(img(:))));
stim_img_perfect =reshape(stim_img_perfect,gridSzX,gridSzY);

img = reshape(img(:),gridSzX,gridSzY);

figure;
subplot(1,3,1);
imagesc(img);axis image;colormap gray;caxis([-0.3,0.5]);
hold on;
plot(on.conesY,on.conesX,'r*');
hold on;

for icone=1:on.nCones
    plot(x+on.conesY(icone),(on.conesX(icone)+y),'r');                      % Plot the circle 
    axis equal;
    grid on;
end
hold on;

R = 3*off.coneGaussSd;                         % Unit radius
x = R*cos(angle);  y = R*sin(angle);   % Coordinates of the circle

plot(off.conesY,off.conesX,'b*');
hold on;
n=100;
for icone=1:off.nCones
    plot(x+off.conesY(icone),y+off.conesX(icone),'b');                      % Plot the circle 
    axis equal;
    grid on;
end

subplot(1,3,2);
imagesc(stim_img);
axis image;colormap gray;caxis([-0.3,0.5]);


subplot(1,3,3);
imagesc(stim_img_perfect);
axis image;colormap gray;caxis([-0.3,0.5]);
