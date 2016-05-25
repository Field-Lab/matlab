clear
datarun.names.rrs_neurons_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data000-from-data000_data003_data006_data007_data008_data009/data000-from-data000_data003_data006_data007_data008_data009.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data000-from-data000_data003_data006_data007_data008_data009/data000-from-data000_data003_data006_data007_data008_data009.params';
datarun.names.rrs_sta_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data000-from-data000_data003_data006_data007_data008_data009/data000-from-data000_data003_data006_data007_data008_data009.sta';
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_sta', 1, 'load_all',1);

datarun=load_data(datarun, opt);
indicies = get_cell_indices(datarun, 338);
sta = datarun.stas.stas{indicies};
 
scale = 2; 

[sig_stixels] = significant_stixels(sta, 'select', 'thresh', 'thresh', 2);
figure
temp_rf = rf_from_sta(sta, 'sig_stixels', sig_stixels);
image = norm_image(temp_rf);

for i = 1:size(image,3)
   upsamp_image(:,:,i) =  imresize(image(:,:,i), scale, 'nearest')
end

sta_img = imagesc(upsamp_image);
colormap gray

hold on

the_fit = datarun.stas.fits{indicies};
ctr = the_fit.mean;
rad = the_fit.sd;
[X,Y] = drawEllipse([ctr*scale rad*scale the_fit.angle]);
plot(X,Y,'Color','k', 'linewidth',2 );

hold on
set(gca, 'ydir', 'reverse')

axis image
drawnow
bounds = [15 31 7 19]*scale;
axis(bounds)

axis off

pix_field_width = datarun.stimulus.field_width*datarun.stimulus.stixel_width;
scale_bar_width = 15; %pixels
stix_scale_bar = (scale_bar_width/pix_field_width)*datarun.stimulus.field_width*scale;
plot([bounds(1)+1; bounds(1)+stix_scale_bar+1], [bounds(4)-1; bounds(4)-1], '-k', 'LineWidth', 2)

title({'ON Type 2'; 'Rod-Sensitive'})
%%

clear
datarun.names.rrs_neurons_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data009-from-data000_data003_data006_data007_data008_data009/data009-from-data000_data003_data006_data007_data008_data009.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data009-from-data000_data003_data006_data007_data008_data009/data009-from-data000_data003_data006_data007_data008_data009.params';
datarun.names.rrs_sta_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data009-from-data000_data003_data006_data007_data008_data009/data009-from-data000_data003_data006_data007_data008_data009.sta';
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_sta', 1, 'load_all',1);

datarun=load_data(datarun, opt);
indicies = get_cell_indices(datarun, 338);
sta = datarun.stas.stas{indicies};

[sig_stixels] = significant_stixels(sta, 'select', 'thresh', 'thresh', 2);
figure
temp_rf = rf_from_sta(sta, 'sig_stixels', sig_stixels);
image = norm_image(temp_rf);
scale =2;
for i = 1:size(image,3)
   upsamp_image(:,:,i) =  imresize(image(:,:,i), scale, 'nearest');
end
% maxi = max(max(upsamp_image(:,:,2)));
% mini = min(min(upsamp_image(:,:,2)));
% sta_img = imagesc((upsamp_image(:,:,2)- mini)./(maxi-mini));
sta_img = imagesc(upsamp_image(:,:,2));

colormap gray


hold on

% the_fit = datarun.stas.fits{indicies};
% ctr = the_fit.mean;
% rad = the_fit.sd;
% [X,Y] = drawEllipse([ctr rad the_fit.angle]);
% plot(X,Y,'Color','k');
% hold on
set(gca, 'ydir', 'reverse')

axis image
drawnow
bounds = [15 31 7 19]*scale;
axis(bounds)

axis off

pix_field_width = datarun.stimulus.field_width*datarun.stimulus.stixel_width;
scale_bar_width = 15; %pixels
stix_scale_bar = (scale_bar_width/pix_field_width)*datarun.stimulus.field_width*scale;
plot([bounds(1)+1; bounds(1)+stix_scale_bar+1], [bounds(4)-1; bounds(4)-1], '-k', 'LineWidth', 2)

title({'ON Type 2'; 'Cone-Sensitive'})
m = 0.25;
caxis([m-2 m+2])

%%

clear
datarun.names.rrs_neurons_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data000-from-data000_data003_data006_data007_data008_data009/data000-from-data000_data003_data006_data007_data008_data009.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data000-from-data000_data003_data006_data007_data008_data009/data000-from-data000_data003_data006_data007_data008_data009.params';
datarun.names.rrs_sta_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data000-from-data000_data003_data006_data007_data008_data009/data000-from-data000_data003_data006_data007_data008_data009.sta';
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_sta', 1, 'load_all',1);

datarun=load_data(datarun, opt);
indicies = get_cell_indices(datarun, 607);
sta = datarun.stas.stas{indicies};
scale =2;

[sig_stixels] = significant_stixels(sta, 'select', 'thresh', 'thresh', 2);
figure
temp_rf = rf_from_sta(sta, 'sig_stixels', sig_stixels);
image = norm_image(temp_rf);

for i = 1:size(image,3)
   upsamp_image(:,:,i) =  imresize(image(:,:,i), scale, 'nearest');
end

sta_img = imagesc(upsamp_image);

colormap gray
hold on

the_fit = datarun.stas.fits{indicies};
ctr = the_fit.mean;
rad = the_fit.sd;
[X,Y] = drawEllipse([ctr*scale rad*scale the_fit.angle]);
plot(X,Y,'Color','k', 'linewidth',2 );

hold on
set(gca, 'ydir', 'reverse')


axis image
drawnow
bounds = [15 26 5 15]*scale;
axis(bounds)

axis off

pix_field_width = datarun.stimulus.field_width*datarun.stimulus.stixel_width;
scale_bar_width = 15; %pixels
stix_scale_bar = (scale_bar_width/pix_field_width)*datarun.stimulus.field_width*scale;
plot([bounds(1)+1; bounds(1)+stix_scale_bar+1], [bounds(4)-1; bounds(4)-1], '-k', 'LineWidth', 2)

title({'ON Type 1'; 'Rod-Sensitive'})


%%

clear
datarun.names.rrs_neurons_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data009-from-data000_data003_data006_data007_data008_data009/data009-from-data000_data003_data006_data007_data008_data009.neurons';
datarun.names.rrs_params_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data009-from-data000_data003_data006_data007_data008_data009/data009-from-data000_data003_data006_data007_data008_data009.params';
datarun.names.rrs_sta_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data009-from-data000_data003_data006_data007_data008_data009/data009-from-data000_data003_data006_data007_data008_data009.sta';
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_sta', 1, 'load_all',1);

datarun=load_data(datarun, opt);
indicies = get_cell_indices(datarun, 607);
sta = datarun.stas.stas{indicies};
scale =2;

[sig_stixels] = significant_stixels(sta, 'select', 'thresh', 'thresh', 2);
figure
temp_rf = rf_from_sta(sta, 'sig_stixels', sig_stixels);
% image = norm_image(temp_rf);
image = temp_rf;

for i = 1:size(image,3)
   upsamp_image(:,:,i) =  imresize(image(:,:,i), scale, 'nearest');
end
maxi = max(max(upsamp_image(:,:,2)));
mini = min(min(upsamp_image(:,:,2)));
sta_img = imagesc((upsamp_image(:,:,2)-mini)./(maxi-mini));
colormap gray
caxis([-1 1])



hold on

the_fit = datarun.stas.fits{indicies};
ctr = the_fit.mean;
rad = the_fit.sd;
[X,Y] = drawEllipse([ctr*scale rad*scale the_fit.angle]);
plot(X,Y,'Color','k', 'linewidth',2 );

hold on
set(gca, 'ydir', 'reverse')


axis image
drawnow
bounds = [15 26 5 15]*scale;
axis(bounds)

axis off

pix_field_width = datarun.stimulus.field_width*datarun.stimulus.stixel_width;
scale_bar_width = 15; %pixels
stix_scale_bar = (scale_bar_width/pix_field_width)*datarun.stimulus.field_width*scale;
plot([bounds(1)+1; bounds(1)+stix_scale_bar+1], [bounds(4)-1; bounds(4)-1], '-k', 'LineWidth', 2)

    title({'ON Type 1'; 'Cone-Sensitive'})


