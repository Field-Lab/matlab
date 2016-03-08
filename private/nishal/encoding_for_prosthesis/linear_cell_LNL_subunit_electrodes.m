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
R =1*on.coneGaussSd;                         % Unit radius
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

R = 1*off.coneGaussSd;                         % Unit radius
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

elecSpacing=2%6;
arrSz=round(spatial_extent*2/elecSpacing);
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
su_elec = full(elecs.su_elec);

nl_ec = @(x) 1./(1+exp(-(x)));
nl_ec_deri = @(x) (exp(-(x))./(1+exp(-(x))).^2);
figure;plot([-3:0.1:7],nl_ec([-3:0.1:7]));hold on;
plot([-3:0.1:7],nl_ec_deri([-3:0.1:7]));

%cell_current = @(current) weight_elec_su*nl_ec(su_elec*current);
%cell_current_linear = @(current, current_old) (cell_current(current_old) + weight_elec_su*(nl_ec_deri(su_elec*current_old).*(su_elec*(current-current_old))));

%cell_current = @(current) nl_ec(weight_elec_su*nl_ec(su_elec*current));
%cell_current_linear = @(current, current_old) (cell_current(current_old) ...
%    + nl_ec_deri(weight_elec_su*(nl_ec(su_elec*current_old))).*  (weight_elec_su*(nl_ec_deri(su_elec*current_old).*(su_elec*(current-current_old)))));

nl1 = @(x) (log(1+exp(x)));
nl2 = @(x) (nl_ec(9*x-4.5));
cell_current = @(current) nl2(weight_elec_su*nl1(su_elec*current));

%% To be modified after this --

%% save architecture for pythonizing! 

save('/Volumes/Lab/Users/bhaishahster/encoding_for_prosthesis_simulation.mat','weight_elec_su','su_elec','stas','stas_inv','gridSzX','gridSzY','-v7.3');

%% %% Input image
img = imread('/Volumes/Lab/Users/bhaishahster/SIPI database/misc/4.2.03.tiff'); 
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
electrodes_chosen = [263,264];%iidx(elec_list>0);
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
current_lim = [-20:2:20];
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

%% Plot electrodes, cells , sub-units

figure;plot(elecs.x,elecs.y,'r.','MarkerSize',10);hold on;plot(elecs.suX,elecs.suY,'b.');
cols = distinguishable_colors(on.nCones);
for icell=1:on.nCones
   iidx_su = 1:elecs.nSU;
    su_list = iidx_su(logical(on.elecs.weight_elec_su(icell,:)));
   for isu=su_list
       hold on;
       plot([on.conesX(icell),elecs.suX(isu)],[on.conesY(icell),elecs.suY(isu)],'Color',cols(icell,:));
       
       
   end
end
hold on;
plot(on.conesX(35),on.conesY(35),'.','MarkerSize',40);
% plot SU
for isu=1:elecs.nSU
    iidx_elecs = 1:elecs.nCones;
    elecs_list = iidx_elecs(logical(elecs.su_elec(isu,:)~=0));
    a = convhull(elecs.x(elecs_list),elecs.y(elecs_list));
    hold on;
    plot(elecs.x(elecs_list(a)),elecs.y(elecs_list(a)),'g');
end

for ielec=1:elecs.nCones
    hold on;
text(elecs.x(ielec),elecs.y(ielec),sprintf('%d',ielec));
end

%% plot electrodes, cells
figure;
%plot(elecs.x,elecs.y,'k.','MarkerSize',10);
%xlim([0,64]);ylim([0,64]);
set(gca,'xTick',[]);set(gca,'yTick',[]);
set(gca,'visible','off');
axis square 

n=100;
angle = 0:2*pi/n:2*pi;            % vector of angles at which points are drawn
R =1*on.coneGaussSd;                         % Unit radius
x = R*cos(angle);  y = R*sin(angle);   % Coordinates of the circle

%plot(on.conesX,on.conesY,'r*');
hold on;

for icone=1:on.nCones
    plot(x+on.conesY(icone),y+64-on.conesX(icone),'r','LineWidth',2);                      % Plot the circle 
    axis square;
    grid on;
end
hold on;
xlim([-1,66]);ylim([-1,66]);



figure;
%plot(elecs.x,elecs.y,'k.','MarkerSize',10);
%xlim([0,64]);ylim([0,64]);
set(gca,'xTick',[]);set(gca,'yTick',[]);
set(gca,'visible','off');
axis square 

n=100;
angle = 0:2*pi/n:2*pi;            % vector of angles at which points are drawn
R =1*off.coneGaussSd;                         % Unit radius
x = R*cos(angle);  y = R*sin(angle);   % Coordinates of the circle

%plot(on.conesX,on.conesY,'r*');
hold on;

for icone=1:on.nCones
    plot(x+off.conesY(icone),y+64-off.conesX(icone),'b','LineWidth',2);                      % Plot the circle 
    axis square;
    grid on;
end
hold on;
xlim([-1,66]);ylim([-1,66]);

% ON and OFF grid without the electrodes!

%% Load theano results
dat = load('/Volumes/Lab/Users/bhaishahster/encoding_for_prosthesis_simulation_perception.mat');
figure;
imagesc(dat.img_recons2);
axis image
colormap gray;
set(gca,'visible','off')

figure;
imagesc(dat.img_perfect2);
axis image
colormap gray;
set(gca,'visible','off')


figure;
imagesc(dat.img_orig);
axis image
colormap gray;
set(gca,'visible','off')

figure;
current = dat.current;
scatter(elecs.y,64-elecs.x,100,20*(current),'filled');colormap gray
%hold on;
%scatter(elecs.y,64-elecs.x,40*abs((current)'+1),zeros(nElecs,1));
axis square;
hold on;
xlim([1,64]);ylim([1,64]);
set(gca,'visible','off')


% second sight current input
figure;
current = current_inp;
scatter(elecs.y,64-elecs.x,100,30*(current),'filled');colormap gray
%hold on;
%scatter(elecs.y,64-elecs.x,40*abs((current)'+1),zeros(nElecs,1));
axis square;
hold on;
xlim([1,64]);ylim([1,64]);
set(gca,'visible','off')