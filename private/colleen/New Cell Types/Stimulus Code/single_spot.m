

clear
%% ------------------------------ INPUTS -----------------------------------
cells = {876};
file_name = '2006-06-06-0/data000-mg/data000/data000';
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-16-1-0.48-11111.xml';
file_path = ['/Users/colleen/matlab/private/colleen/New Cell Types/Stimulus Code/2006-06-06-0/data000/'];
screen_width = 640; % in pixels % vertical size
screen_height = 320; % hortizontal size
stixels_ref = 16;
stixels_focal = 1; % must be a factor of stixels_ref
%% ------------------------------- Load Data ------------------------------------------

if ~exist(file_path)
    mkdir(file_path)
end

datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];
opt=struct('verbose',1,'load_params',1,'load_neurons',0,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',false);
opt.load_sta_params.frames = 1:30;% if this line is missing, will error; have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

%% ------------------------------- Plot Vision STA -----------------------------
myMap = zeros(screen_height, screen_width); % pixesl on the screen
[cell_numbers, cell_type, cell_type_number] = get_cell_indices(datarun, cells{1});
sta = datarun.stas.stas{cell_numbers};
%     sig_stixels = significant_stixels(sta);
sig_stixels = significant_stixels(sta)

% find peak frame
time_course = time_course_from_sta(sta, sig_stixels);


temp_rf = rf_from_sta(sta, 'sig_stixels', sig_stixels);

the_fit = datarun.stas.fits{cell_numbers};
ctr = the_fit.mean;
rad = the_fit.sd;
[X,Y] = drawEllipse_upsampled([ctr rad the_fit.angle]);
%         Y = screen_height/stixels_ref - Y;
%         figure;
%                 plot(X,Y,'Color','k');
axis([0 screen_width/stixels_ref 0 screen_height/stixels_ref])

X_large =  round(screen_width*X/(screen_width/stixels_ref));
Y_large =  round(screen_height*Y/(screen_height/stixels_ref));
%        ellipse_t = fit_ellipse( X_large,Y_large)
%        myMap([Y_large; X_large]') = 1;
for i = 1:length(X_large)
    myMap(Y_large(i),X_large(i)) = 1;
end

myMap_filled = imfill(myMap,'holes');


figure
imagesc(myMap_filled)
axis equal

%% ------------------------------- Write the mask to a file ---------------------------------
dlmwrite([file_path,num2str(cells{1}), '.txt'], myMap_filled, 'delimiter', '\t', 'newline', 'pc'); % if this errors which save path
savedMap = dlmread(file_path);

%% ------------------------------- Display mask ----------------------------------------------
figure
imagesc(savedMap)
title('Stixels to be modulated')
axis equal