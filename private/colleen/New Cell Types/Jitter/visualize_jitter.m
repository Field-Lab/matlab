close all
clear

date = '2015-09-23-7';
concate = 'data031';
% cell_ids = [81 530 532 918 1006 1010 1103 1175 1597 1880 2378 2421 2644 3394 3398 3559 3931 4043 4132 4341 4658 4850 4864 5075 5090 5299 5342 5361 5493 5506 6323 6788 6830 7087 7131 7145 7203];

cell_ids =  [185 306 413 996 1058 1759 1821 2061 3396 3830 4656 5342 5733 6064 6213 6831 7148 7628 2613 3215 532 1022 1568 2376 3186 4041 4703 5239 6788 7101 917 1026 3931 7216]
savepath = ['/Volumes/Lab/Users/crhoades/Jitter/', date, '/', concate];
if ~exist(savepath)
    mkdir(savepath);
end

time_course_total = [];
summed_sta =zeros(320,640,3);
for n= 1:length(cell_ids)
    clearvars sta temp
    filename = ['/Volumes/Lab/Users/crhoades/Jitter/', date, '/',concate, '/Cell ',num2str(cell_ids(n)), '.mat'];
    load(filename);
    %     if ~exist('sta')
    try
     sta = temp;
    catch
    end
    
    %     end
    
    
    %loads sta variable
    
    sig_stixels = significant_stixels(sta, 'select', 'thresh', 'thresh', 3);
    
    width = size(sta, 2);
    height = size(sta,1);
    % full = full(sig_stixels);
    border_width = 10;
    sig_stixels([1:border_width, (height - border_width+1):height], :) = 0;
    sig_stixels(:,[1:border_width, (width - border_width+1):width]) = 0;
    
    % sig_stixels_new = sparse(full);
    [timecourse, params] = time_course_from_sta(sta, sig_stixels);
    
    figure;
    subplot(2,1,2)
    set(gca, 'ColorOrder', eye(3), 'NextPlot', 'replacechildren');
    co = get(gca, 'ColorOrder');
    plot(repmat(linspace(-44*8.333,0,45)',1,3),timecourse)
    axis tight
    
    timecourse = timecourse./max(abs(timecourse(:)));
    time_course_total = [time_course_total, timecourse];
    
    [a,b] = max(abs(timecourse));
    peak_frame = max(b);
    x=sta;
    mx=max(x(:));
    mn=min(x(:));
    upper = 65535; %16 bit
    lower  =0;
    y = (x - mn)*(upper-lower)/(mx - mn) +lower;
    y16 = uint16(y - 1);
    subplot(2,1,1)
    % for i = 1:size(sta,4)
    %     imagesc(y16(:,:,:,i));
    %     pause(0.5)
    % end
    imagesc(y16(:,:,:,peak_frame));
    
    summed_sta = summed_sta + sta(:,:,:,peak_frame);
    %
    % rf = rf_from_sta(sta, 'sig_stixels', sig_stixels);
    % figure;
    % if sum(timecourse(b)) < 0
    %     imagesc(norm_image(-rf));
    % else
    %         imagesc(norm_image(rf));
    % end
    axis equal
    axis off
    
    title({date; concate; ['Cell ', num2str(cell_ids(n))]})
    savefig(gcf, [savepath, '/Cell ' num2str(cell_ids(n))])
end

figure;
set(gca, 'ColorOrder', repmat(eye(3), length(cell_ids),1), 'NextPlot', 'replacechildren');
co = get(gca, 'ColorOrder');
plot(repmat(linspace(-44*8.333,0,45)',1, length(cell_ids)*3),time_course_total,'--')
axis tight
title({date; concate; ['TimeCourses']})
hold on
% set(gca, 'ColorOrder', repmat(eye(3), length(cell_ids),1), 'NextPlot', 'replacechildren');
% co = get(gca, 'ColorOrder');
plot(repmat(linspace(-44*8.333,0,45)',1, (length(cell_ids)-9)*3),time_course_total(:,28:end))
axis tight

savefig(gcf, [savepath, 'Timecourses LBC SBC'])

%     x=sta;
%     mx=max(summed_sta(:));
%     mn=min(summed_sta(:));
%     upper = 65535; %16 bit
%     lower  =0;
%     y = (summed_sta - mn)*(upper-lower)/(mx - mn) +lower;
%     y16 = uint16(y - 1);
%
% figure; imagesc(y16);



%% compare to Vision STA
% clear
% close all
dbstop if error
dataparam.date='2015-09-23-7';
dataparam.concatname='d19-39/data030-from-data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-20-8-0.48-11111.xml';
dataparam.file_name_right = [dataparam.date, '/', dataparam.concatname,'/data030-from-data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039'];


% Right Movie
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name_right, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name_right, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name_right, '.sta'];
savepath = ['/Volumes/Lab/Users/crhoades/Jitter/RegularSTAs/', dataparam.date, '/', dataparam.concatname];
if ~exist(savepath)
    mkdir(savepath);
end

% cell_ids = [81 530 532 918 1006 1010 1103 1175 1597 1880 2378 2421 2644 3394 3398 3559 3931 4043 4132 4341 4658 4850 4864 5075 5090 5299 5342 5361 5493 5506 6323 6788 6830 7087 7131 7145 7203];
% 
% Load Data2
slashes = strfind(datarun.names.rrs_neurons_path, '/');
dataset = datarun.names.rrs_neurons_path(slashes(3)+1:slashes(5)-1);
to_replace = strfind(dataset, '/');
dataparam.dataset(to_replace) = '-';

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);
time_course_total = [];
cell_indices = get_cell_indices(datarun,cell_ids);
for i = 1:length(cell_ids)
    sta = datarun.stas.stas{cell_indices(i)};
    
    % sig_stixels_new = sparse(full);
    
    sig_stixels = significant_stixels(sta, 'select', 'thresh', 'thresh', 3);
    [timecourse, params] = time_course_from_sta(sta, sig_stixels);
    
    %     timecourse = datarun.stas.time_courses{cell_ids(i)};
    
    figure;
    subplot(2,1,2)
    set(gca, 'ColorOrder', eye(3), 'NextPlot', 'replacechildren');
    co = get(gca, 'ColorOrder');
    plot(repmat(linspace(-(length(timecourse)- 1)*8.333*dataparam.mdf_file(end-10),0,length(timecourse))',1,3),timecourse)
    axis tight
    
    timecourse = timecourse./max(abs(timecourse(:)));
    time_course_total = [time_course_total, timecourse];
    
    [a,b] = max(abs(timecourse));
    peak_frame = max(b);
    x=sta;
    mx=max(x(:));
    mn=min(x(:));
    upper = 65535; %16 bit
    lower  =0;
    y = (x - mn)*(upper-lower)/(mx - mn) +lower;
    y16 = uint16(y - 1);
    subplot(2,1,1)
    % for i = 1:size(sta,4)
    %     imagesc(y16(:,:,:,i));
    %     pause(0.5)
    % end
    imagesc(y16(:,:,:,peak_frame));
    axis equal
    axis off
    title({dataset; ['Cell ',num2str(cell_ids(i))]});
    savefig(gcf, [savepath, '/Cell ' num2str(cell_ids(i))])
end

figure;
set(gca, 'ColorOrder', repmat(eye(3), length(cell_ids),1), 'NextPlot', 'replacechildren');
co = get(gca, 'ColorOrder');
plot(repmat(linspace(-(size(timecourse,1)-1)*8.333*dataparam.mdf_file(end-10),0,size(timecourse,1))',1, length(cell_ids)*3),time_course_total)
axis tight
title({dataparam.date; dataparam.concatname(1:14); ['TimeCourses']})


