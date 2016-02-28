%% Linear cells, LNL sub-unit stimulation model

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
R = 2*on.coneGaussSd;                         % Unit radius
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

R = 2*off.coneGaussSd;                         % Unit radius
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

elecSpacing=4%6;
arrSz=round(spatial_extent*1.7/elecSpacing);
elecLatticeOrientation = pi/6;
elecs = getElectrodes_simulation(elecSpacing,arrSz,elecLatticeOrientation,spatial_extent)
nElecs = length(elecs.x);
    
figure;plot(elecs.x,elecs.y,'r.','MarkerSize',10);hold on;plot(elecs.suX,elecs.suY,'b.');
on = map_population_electrodes_subunits(on,elecs)
off = map_population_electrodes_subunits(off,elecs)

%% Experimentation

gridSzX =  on.gridSzX;
gridSzY = on.gridSzY;

stas_on = on.stas; 
stas_off = off.stas;
stas = [stas_on;-stas_off];
stas_inv = pinv(stas);
nCells = size(stas,1);


weight_elec_su = [on.elecs.weight_elec_su;off.elecs.weight_elec_su];
su_elec = elecs.su_elec;

nl_ec = @(x) 1./(1+exp(-(x)));
nl_ec_deri = @(x) (exp(-(x))./(1+exp(-(x))).^2);
figure;plot([-3:0.1:7],nl_ec([-3:0.1:7]));hold on;
plot([-3:0.1:7],nl_ec_deri([-3:0.1:7]));

%cell_current = @(current) weight_elec_su*nl_ec(su_elec*current);
%cell_current_linear = @(current, current_old) (cell_current(current_old) + weight_elec_su*(nl_ec_deri(su_elec*current_old).*(su_elec*(current-current_old))));

cell_current = @(current) nl_ec(weight_elec_su*nl_ec(su_elec*current));
cell_current_linear = @(current, current_old) (cell_current(current_old) ...
    + nl_ec_deri(weight_elec_su*(nl_ec(su_elec*current_old))).*  (weight_elec_su*(nl_ec_deri(su_elec*current_old).*(su_elec*(current-current_old)))));

%% To be modified after this --
%% %% Input image
img = imread('~/Downloads/SIPI database/misc/4.2.03.tiff'); 
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
obj_log=[];obj_log2 = [];

figure;
for iter=1:50
    iter
cvx_begin quiet
variables cell_r_norm(nCells) current(nElecs)
minimize ((sum_square(stas_inv*(cell_r_norm*b+a) - img_flat)) + 1*sum_square(current-current_old))

subject to 
 cell_current_linear(current,current_old) == cell_r_norm
   cell_r_norm<=1
   cell_r_norm>=0
cvx_end
current_old=current;
obj_log=[obj_log;((sum_square(stas_inv*(cell_r_norm*b+a) - img_flat)))];
obj_log2 = [obj_log2;(sum_square(stas_inv*(cell_r_norm*b+a) - img_flat)) + 0.1*sum_square(current-current_old)];
hold on;plot(current);
end

figure;
plot(obj_log);
figure;
plot(obj_log2);
cell_r_norm = cell_current(current);
stim_img = stas_inv*(cell_r_norm*b+a);
stim_img =reshape(stim_img,gridSzX,gridSzY);

figure;
hold on;
scatter(elecs.x,elecs.y,40*abs((current)+1),zeros(nElecs,1));
axis equal
hold on;
imagesc(reshape(img_flat,gridSzX,gridSzY));colormap gray
for icell=1:on.nCones
    hold on;
    pos = [on.conesX(icell)-2*on.coneGaussSd,on.conesY(icell)-2*on.coneGaussSd,4*on.coneGaussSd,4*on.coneGaussSd];
    rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[1,1,1]*cell_r_norm(icell))
    axis equal
end
%hold on;
%plot(on.conesX,on.conesY,'r*');
hold on;
scatter(elecs.x,elecs.y,40*abs((current)+1),2*sign(current),'filled');colormap cool
hold on;
scatter(elecs.x,elecs.y,40*abs((current)+1),zeros(nElecs,1));
axis equal


figure;
%hold on;
%scatter(elecs.x,elecs.y,40*abs((current)+1),zeros(nElecs,1));
%axis equal
%hold on;
%imagesc(reshape(img_flat,gridSzX,gridSzY));colormap gray
for icell=1:off.nCones
    hold on;
    pos = [off.conesX(icell)-2*off.coneGaussSd,off.conesY(icell)-2*off.coneGaussSd,4*off.coneGaussSd,4*off.coneGaussSd];
    rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[1,1,1]*cell_r_norm(icell+on.nCones))
    axis equal
end
%hold on;
%plot(on.conesX,on.conesY,'r*');
hold on;
scatter(elecs.x,elecs.y,40*abs((current)+1),2*sign(current),'filled');colormap cool
hold on;
scatter(elecs.x,elecs.y,40*abs((current)+1),zeros(nElecs,1));
axis equal


%% Show reconstruction
figure;

subplot(1,5,1);
imagesc(reshape(img_flat,gridSzX,gridSzY));axis image;colormap gray
title('Input image');

subplot(1,5,2);
imagesc(stim_img);axis image;
title('Achieved Stimulation');
colormap gray

% Perfect stimulation
subplot(1,5,3);
perfect_cell_resp = stas*(img_flat);
stim_img_perfect = stas_inv*(perfect_cell_resp);
stim_img_perfect =reshape(stim_img_perfect,gridSzX,gridSzY);
imagesc(stim_img_perfect);axis image
colormap gray;
title('Perfect stimulation');

% second sight stimulation 
elec_sta_inp = elecs.stas*(img_flat);
current_inp = elec_sta_inp; % ??
cell_r_norm_pros =cell_current(current_inp);
stim_img_pros = stas_inv*(cell_r_norm_pros*b+a);
stim_img_pros =reshape(stim_img_pros,gridSzX,gridSzY);
subplot(1,5,4);
imagesc(stim_img_pros);axis image
colormap gray;
title('Current Prosthesis');
% 
% % second sight stimulation with only ON cells
% elec_sta_inp = elecs.stas*(img_flat);
% current_inp = elec_sta_inp; % ??
% cell_r_norm_pros =cell_current(current_inp);
% %cell_r_norm_pros_off = cell_current(0*current_inp);
% cell_r_norm_pros(on.nCones+1:on.nCones+off.nCones)=3;
% stim_img_pros = stas_inv*(cell_r_norm_pros_off*b+a);
% stim_img_pros =reshape(stim_img_pros,gridSzX,gridSzY);
% subplot(1,5,5);
% imagesc(stim_img_pros);axis image
% colormap gray;
% title('Current Prosthesis');

%% bi-electrode stimulation . piece wise linear curves?

cellProbe = zeros(on.nCones+off.nCones,1);
cell_chosen=35;
cellProbe(cell_chosen)=1;

weight_elec_su = [on.elecs.weight_elec_su;off.elecs.weight_elec_su];
su_elec = elecs.su_elec;
elec_list = su_elec'*( weight_elec_su'*cellProbe);
iidx = 1:nElecs;
electrodes_chosen = iidx(elec_list>0);
figure;
nElec_chosen = length(electrodes_chosen);
iiicnt=0;


for iielec=1:nElec_chosen
    for jjelec=1:nElec_chosen
        iiicnt=iiicnt+1;
        
        subplot(nElec_chosen,nElec_chosen,iiicnt);
        ielec = electrodes_chosen(iielec);
        jelec = electrodes_chosen(jjelec);
        if(ielec==jelec)
        continue;    
        end

icnt=0;
resp = [];
current_lim = [-10:0.5:10];
for icurrent = current_lim;
    icnt=icnt+1;
    jcnt=0;
    for jcurrent = current_lim;
        jcnt=jcnt+1;
        current_inp=-10*ones(nElecs,1);
        current_inp(ielec)=icurrent;
        current_inp(jelec)=jcurrent;
        
        cell_r_norm_bielec =cell_current(current_inp);
        resp(icnt,jcnt)=cell_r_norm_bielec(cell_chosen);
    end
end

nrep = length(current_lim);

contourf(repelem(current_lim,nrep,1),repelem(current_lim',1,nrep),resp,30);
%imagesc(resp<0.5)
hold on;
plot(zeros(nrep,1),current_lim,'k');
hold on;
plot(current_lim,zeros(nrep,1),'k');
set(gca,'xTick',[]);
set(gca,'yTick',[]);
    end
end