elec_sta_inp = elecs.stas*(img_flat);
current_inp = elec_sta_inp; % ??
cell_r_norm_pros =cell_current(current_inp);
stim_img_pros = stas_inv*(cell_r_norm_pros*b+a);
stim_img_pros =reshape(stim_img_pros,gridSzX,gridSzY);
subplot(1,5,4);
imagesc(stim_img_pros);axis image
colormap gray;
title('Current Prosthesis');

%% Second sight prosthesis
cell_current = @(x)1./(1+exp(-(1.46*((weight_elec_su*(su_elec*x)))-2.19)));
elec_sta_inp = elecs.stas*(img_flat)/4;
current_inp = elec_sta_inp; % ??
cell_r_norm_pros =cell_current(current_inp);
%cell_r_norm_pros(1:on.nCones)= cell_r_norm_pros(1:on.nCones)*1;
%cell_r_norm_pros(on.nCones+1:end)= cell_r_norm_pros(on.nCones+1:end)*0;
%stim_img_pros = stas_inv*(cell_r_norm_pros*b+a);
%stim_img_pros = pinv(on.stas)*(cell_r_norm_pros(1:on.nCones)*b+a);
scale_on=1.2;scale_off=0.8; % 1.2, or 0.8
stim_img_pros = (on.stas') *(cell_r_norm_pros(1:on.nCones)*b+a)*scale_on - scale_off*(off.stas') *(cell_r_norm_pros(on.nCones+1:end)*b+a);
stim_img_pros =reshape(stim_img_pros,gridSzX,gridSzY);
figure;
hist(cell_r_norm_pros);
figure;
imagesc(stim_img_pros);axis image
colormap gray;
set(gca,'visible','off')
title('Current Prosthesis');

figure;
current = current_inp;
scatter(elecs.y,64-elecs.x,100,30*(current),'filled');colormap gray
%hold on;
%scatter(elecs.y,64-elecs.x,40*abs((current)'+1),zeros(nElecs,1));
axis square;
hold on;
xlim([1,64]);ylim([1,64]);
set(gca,'visible','off')



figure;
%hold on;
%scatter(elecs.x,elecs.y,40*abs((current)+1),zeros(nElecs,1));
%axis equal
%hold on;
%imagesc(reshape(img_flat,gridSzX,gridSzY));colormap gray
for icell=1:off.nCones
    hold on;
    pos = [off.conesX(icell)-2*off.coneGaussSd,off.conesY(icell)-2*off.coneGaussSd,2*off.coneGaussSd,2*off.coneGaussSd];
    rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[1,1,1]*cell_r_norm_pros(icell+on.nCones))
    axis equal
end
% %hold on;
% %plot(on.conesX,on.conesY,'r*');
% hold on;
% scatter(elecs.x,elecs.y,40*abs((current)+1),2*sign(current),'filled');colormap cool
% hold on;
% scatter(elecs.x,elecs.y,40*abs((current)+1),zeros(nElecs,1));
% axis equal

figure;
% hold on;
% scatter(elecs.x,elecs.y,40*abs((current)+1),zeros(nElecs,1));
axis equal
hold on;
imagesc(reshape(img_flat,gridSzX,gridSzY));colormap gray
for icell=1:on.nCones
    hold on;
    pos = [on.conesX(icell)-2*on.coneGaussSd,on.conesY(icell)-2*on.coneGaussSd,2*on.coneGaussSd,2*on.coneGaussSd];
    rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[1,1,1]*cell_r_norm_pros(icell))
    axis equal
end
figure;

imagesc(reshape(img_flat,gridSzX,gridSzY));colormap gray
axis image