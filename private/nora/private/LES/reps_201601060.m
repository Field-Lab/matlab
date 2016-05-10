%% Dataruns: run first for everything

clear
regA{1} = 'data003'; 
sigma2{1} = 'data006';
regB{1} = 'data007'; 
sigma4{1} = 'data008';
cells{1} = 'Maskin2';

regA{2} = 'data007'; %cells = [275;2343;3483;4420;5717];
sigma2{2} = 'data010'; %cells = [345;2464;3379;4413;5716];
regB{2} = 'data011'; %cells = [275;2341;3377;4412;5719];
sigma4{2} = 'data012'; %cells = [202;2465;3481;4411;5717];
cells{2} = 'Maskin';

classification = 'data001';

Analysis_Path = '/Volumes/Analysis/2016-01-05-0/map-from-data001/';
fig_save = '/Users/Nora/Desktop/Fig_Output';
%mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

% classification run
class_datarun = load_data(['/Volumes/Analysis/2016-01-05-0/streamed/' classification '/' classification]);
class_datarun = load_params(class_datarun);
class_datarun = load_neurons(class_datarun);

GS = 0;


%% check stability

for maskin = 1:2;
    
    datarun = load_data([ Analysis_Path regA{maskin} '/' regA{maskin}]);
    datarun = load_neurons(datarun);
    regA_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells{maskin}, 'datarun_class', class_datarun, 'visual_check', 0);
    datarun = load_data([Analysis_Path regB{maskin} '/' regB{maskin}]);
    datarun = load_neurons(datarun);
    regB_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells{maskin}, 'datarun_class', class_datarun, 'visual_check', 0);
    
    for i = 1:4
        try
            figure;
            IDP_plot_PSTH(regA_data, i);
            hold on; IDP_plot_PSTH(regB_data, i);
        catch
        end
    end
    
end

%% load NSinterval: only for GS!
if GS
    i_chunk = 1;
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
    NSmovie = movie_chunk;
    i_chunk = 2;
    
    % mask movie
    while size(NSmovie,3) < 1200
        load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
        NSmovie = cat(3,NSmovie, movie_chunk);
        if size(NSmovie,3) > 1200
            NSmovie = NSmovie(:,:,1:1200);
        end
        i_chunk = i_chunk + 1;
    end
    NSmovie = imresize(NSmovie, 1/8, 'box');
    NSmovie = permute(NSmovie, [2 1 3]);
    NSmovie = NSmovie - 64;
end

%% load LES movie: only for GS
if GS
    maskin = 1; sigma = 2;
    class_datarun = load_sta(class_datarun, 'load_sta', cells{maskin}, 'save_rf', true);
    if sigma == 2
        datarun = load_data([Analysis_Path sigma2{maskin} '/' sigma2{maskin}]);
        if maskin == 1
            mask_file = load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/testmask_1922_stix2/mask_and_filter.mat');
        else
            mask_file = load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/mask_320_sigma2/mask_and_filter.mat');
        end
    else
        datarun = load_data([Analysis_Path sigma4{maskin} '/' sigma4{maskin}]);
        if maskin == 1
            mask_file = load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/testmask_1922_sigma4/mask_and_filter.mat');
        else
            mask_file = load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/mask_320_sigma4/mask_and_filter.mat');
        end
    end
    maskin = 1; sigma = 2;
    if sigma == 2
        if maskin == 1
            movie_file = 'testmask_1922_stix2';
        else
            movie_file = 'mask_320_sigma2';
        end
    else
        if maskin == 1
            movie_file = 'testmask_1922_sigma4';
        else
            movie_file = 'mask_320_sigma4';
        end
    end
    i_chunk = 1;
    load(['/Volumes/Lab/Users/Nora/new_stim_nora/mask_NSEM/masks/' movie_file '/LES/movie_chunk_' num2str(i_chunk) '.mat']);
    LESmovie = movie_chunk;
    i_chunk = 2;
    while size(LESmovie,3) < 1200
        load(['/Volumes/Lab/Users/Nora/new_stim_nora/mask_NSEM/masks/' movie_file '/LES/movie_chunk_' num2str(i_chunk) '.mat']);
        LESmovie = cat(3,LESmovie, movie_chunk);
        if size(LESmovie,3) > 1200
            LESmovie = LESmovie(:,:,1:1200);
        end
        i_chunk = i_chunk + 1;
    end
    
    
    LESmovie = imresize(LESmovie, 1/8, 'box');
    LESmovie = permute(LESmovie, [2 1 3]);
    LESmovie = LESmovie - 64;

    mask = imresize(mask_file.final_mask, 1/8, 'box');
    mask = repmat(mask, 1, 1, 1200);
    mask_movie = NSmovie .* mask;
    LESmovie = LESmovie .* mask;
end

%% look at STA versus filter versus mask size
set(0, 'defaultFigurePosition', [2 2 10 8])
maskin = 1; sigma = 2;
class_datarun = load_sta(class_datarun, 'load_sta', 'all', 'save_rf', true);
ids = get_cell_indices(class_datarun, cells{maskin});
% load masks
if sigma == 2
    datarun = load_data([Analysis_Path sigma2{maskin} '/' sigma2{maskin}]);
    if maskin == 1
        mask_file = load('/Users/Nora/Desktop/Data/2016-01-05-0-Visual/2016-01-05-0/testmask_1922_stix2/mask_and_filter.mat');
    else
        mask_file = load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/mask_320_sigma2/mask_and_filter.mat');
    end
else
    datarun = load_data([Analysis_Path sigma4{maskin} '/' sigma4{maskin}]);
    if maskin == 1
        mask_file = load('/Users/Nora/Desktop/Data/2016-01-05-0-Visual/2016-01-05-0/testmask_1922_sigma4/mask_and_filter.mat');
    else
        mask_file = load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/mask_320_sigma4/mask_and_filter.mat');
    end
end

% i_cell = 1;
% % plot STA with mask border
% sta = imresize(sum(class_datarun.stas.rfs{ids(i_cell)},3), 8, 'nearest');
% figure(1); imagesc(sta + edge(mask_file.final_mask)); axis image
% caxis([min(sta(:)), max(sta(:))])
% colormap gray
% %title([cells{i_cell} ' Cell ' num2str(i_cell) ' 2 Sigma'])
% clear sta
% %figure; imagesc(sum(mask_file.linear_filter,3)+ edge(mask_file.final_mask)); axis image


avg_rf = sum(get_average_rf(class_datarun, 'On Parasol', 'scale', 5),3);
n_angles = 10;
slice = zeros(101,1);
range = -50:50;

for i = 1:4
    %angle = i*360/n_angles;
    %temp = imrotate(avg_rf, angle, 'crop');
    figure(1); imagesc(avg_rf); axis image; title(i)
    prof = improfile;
    prof = [zeros(100,1); prof; zeros(100,1)];
    [~,center] = max(prof);
    %figure(2); hold on; plot(prof(center+range));
    slice = slice + prof(center+range);
end
slice = slice/n_angles;
plot(slice)

a = fit((1:101)'/5,slice, fittype('gauss1'));
width = a.c1/sqrt(2)
avg_profile = [slice, (1:101)'/(5*width)];

for cell_numbers = ids
    the_fit = class_datarun.stas.fits{cell_numbers};
    ctr = the_fit.mean;
    rad = mean(the_fit.sd)
end


%%
% range = -10:10;
% avg_rf = sum(avg_rf,3);
% figure; subplot(1,3,2);
% center = [58 28];
% imagesc(sum(avg_rf(range+center(2),range+center(1),:),3));
% axis image
% for i_cell = [1 3]+1
%     sta = imresize(sum(class_datarun.stas.rfs{ids(i_cell)},3), 3, 'nearest');
%     center = 3*class_datarun.vision.sta_fits{ids(i_cell)}.mean;
%     subplot(1,3,i_cell-1);
%     imagesc(sta(range+(60-round(center(2))),range+round(center(1))));
%     axis image
% end

%% look at the generator signal

if GS
    datarun = load_neurons(datarun);
    sigma_data = interleaved_data_prep(datarun, 1100, 120, 'cell_spec', cells{maskin}, 'datarun_class', class_datarun, 'visual_check', 0);
    
    % separate the conditions
    mask_data.testspikes = sigma_data.testspikes(1:30, :);
    LES_data.testspikes = sigma_data.testspikes(31:60, :);
    %comp_data.testspikes = sigma_data.testspikes(61:90, :);
    %compLES_data.testspikes = sigma_data.testspikes(91:end, :);
    
    % load STAs
    class_datarun = load_sta(class_datarun, 'load_sta', cells{maskin}, 'save_rf', true);
    ids = get_cell_indices(class_datarun, cells{maskin});
    
    for i_cell = 1:4
        
        % plot STA with mask border
        sta = imresize(sum(class_datarun.stas.rfs{ids(i_cell)},3), 8, 'nearest');
        figure(1); imagesc(sta + edge(mask_file.final_mask));
        caxis([min(sta(:)), max(sta(:))])
        colormap gray
        %title([cells{i_cell} ' Cell ' num2str(i_cell) ' 2 Sigma'])
        clear sta
        
        % switch STA for convolution for GS
        if GS
            sta = squeeze(sum(class_datarun.stas.stas{ids(i_cell)},3));
            sta = flip(sta,1);
            sta = flip(sta,2);
            sta = flip(sta,3); 
            figure(2); imagesc(mask_movie(:,:,1));
            figure(3); imagesc(LESmovie(:,:,1));
            figure(4)
            default_colors = get(gca,'ColorOrder');
            response1 = IDP_plot_PSTH(mask_data, i_cell);
            response2 = IDP_plot_PSTH(LES_data, i_cell);
            drive = squeeze(convn(mask_movie, sta, 'valid'));
            %drive = drive+100;
            drive(drive<0) = 0;
            plot(-drive(1:1071), 'Color', default_colors(1,:)-[0 0.2 0.2]);
            hold on
            drive = squeeze(convn(LESmovie,  sta, 'valid'));
            %drive = drive+100;
            drive(drive<0) = 0;
            plot(-drive(1:1071), 'Color', default_colors(2,:)-[0.2 0.2 0]);
            plot(response1(30:1100), 'Color', default_colors(1,:));
            plot(response2(30:1100), 'Color', default_colors(2,:));
            hold off
        end
        
        % exportfig(gcf, [fig_save '/maskin2_sigma2_LESmaskdrive_cell_' num2str(i_cell)], 'Bounds', 'loose', 'Color', 'rgb')
        
        pause()
        %     plot(drive(1:1071), response(30:1100), '.');
        %     [x, y, error, bincount, binedge] = curve_from_binning(drive, response(30:end));
        %     hold on; plot(x,y,'k', 'LineWidth', 2); hold off;
        %     pause()
        %sta_filter = sta_filter + sta.*mask_file{i}.final_mask;
    end
    
end

%imagesc(sta_filter);
%figure; imagesc(sum(mask_file{i}.linear_filter,3));

%% look at all the rasters!

% first 30 are mask
% next 30 are LES
% next 30 are complement
% final 30 are complement + LES
for sigma = 1:2
    % Load the data
    for i = 1:2
        if sigma==1
            datarun = load_data([Analysis_Path sigma2{i} '/' sigma2{i}]);
            fig_offset = 0;
        else
            datarun = load_data([Analysis_Path sigma4{i} '/' sigma4{i}]);
            fig_offset = 3;
        end
        
        if i == 1
            subplot_offset = 0;
        else
            subplot_offset = 4;
        end
        
        datarun = load_neurons(datarun);
        sigma_data = interleaved_data_prep(datarun, 1100, 120, 'cell_spec', cells{i}, 'datarun_class', class_datarun, 'visual_check', 0);
        
        
        % separate the conditions
        mask_data.testspikes = sigma_data.testspikes(1:30, :);
        LES_data.testspikes = sigma_data.testspikes(31:60, :);
        comp_data.testspikes = sigma_data.testspikes(61:90, :);
        compLES_data.testspikes = sigma_data.testspikes(91:end, :);
        
        % reg data
        %         datarun = load_data([ Analysis_Path regA{i} '/' regA{i}]);
        %         datarun = load_neurons(datarun);
        %         regA_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells{i}, 'datarun_class', class_datarun, 'visual_check', 0);
        datarun = load_data([ Analysis_Path regB{i} '/' regB{i}]);
        datarun = load_neurons(datarun);
        regB_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells{i}, 'datarun_class', class_datarun, 'visual_check', 0);
        
        % plot for each cell
        for cell = 1:size(sigma_data.testspikes,2)
            
            % reg versus mask
            %figure(2+fig_offset); subplot(4, 2, cell+subplot_offset); hold on;
            figure(2);
            reg = IDP_plot_PSTH(regB_data, cell);
            mask = IDP_plot_PSTH(mask_data, cell);
            %             xlim([0 1100])
            %             set(gca, 'XTick', []);
            %             set(gca, 'YTick', []);
            %
            %             % reg versus complement
            %             figure(2+fig_offset); subplot(4, 2, cell+subplot_offset); hold on;
            %             plot(reg)
            comp = IDP_plot_PSTH(comp_data, cell);
            %             xlim([0 1100])
            %             set(gca, 'XTick', []);
            %             set(gca, 'YTick', []);
            
            % reg versus complement + mask
            %figure(3+fig_offset); subplot(4, 2, cell+subplot_offset); hold on;
            figure(1);
            plot(reg,'Color', default_colors(4,:))
            hold on
            plot(mask+comp, 'Color', default_colors(5,:))
            plot(-mask, 'Color', default_colors(5,:)+[0.0 0.0 0.3])
            plot(-comp, 'Color', default_colors(5,:)-[0.3 0.3 0.15])
            %xlim([0 1100])
            %set(gca, 'XTick', []);
            %set(gca, 'YTick', []);
            hold off
            exportfig(gcf, [fig_save '/maskin' num2str(mod(i,2)+1) '_sigma' num2str(2*sigma) '_compmask_cell_' num2str(cell)], 'Bounds', 'loose', 'Color', 'rgb')
            
            %             % mask versus LES
            %             figure(4+fig_offset);subplot(4, 2, cell+subplot_offset); hold on;
            %             plot(mask)
            %             LES = IDP_plot_PSTH(LES_data, cell);
            %             xlim([0 1100])
            %             set(gca, 'XTick', []);
            %             set(gca, 'YTick', []);
            
            %             % reg versus LES+complement
            %             figure(5+fig_offset); subplot(4, 2, cell+subplot_offset); hold on;
            %             plot(reg)
            %             plot(LES+comp)
            %             compLES = IDP_plot_PSTH(compLES_data, cell);
            %             xlim([0 1100])
            %             set(gca, 'XTick', []);
            %             set(gca, 'YTick', []);
        end
    end
end

for fig = 1%:10
    figure(fig);
    % exportfig(gcf, [fig_save '/fig' num2str(fig)], 'Bounds', 'loose', 'Color', 'rgb')
end

%% get image transition times
% mask movie
i_chunk = 1;
image_transitions = [1];
while image_transitions(end) < 1200
    load(['/Users/Nora/Desktop/Data/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
    image_transitions = [image_transitions image_transitions(end)+size(movie_chunk,3)];
    i_chunk = i_chunk + 1;
end
image_transitions = image_transitions(1:(end));

%% look at spike sorting

% first 30 are mask
% next 30 are LES
% next 30 are complement
% final 30 are complement + LES
default_colors = get(gca,'ColorOrder');
for sigma = 1%:2
    % Load the data
    for i = 1%:2
        if sigma==1
            datarun = load_data([Analysis_Path sigma2{i} '/' sigma2{i}]);
            datarun_cf = load_data([Analysis_Path sigma2{i} '-cf/' sigma2{i} '-cf']);
            fig_offset = 0;
        else
            datarun = load_data([Analysis_Path sigma4{i} '/' sigma4{i}]);
            datarun_cf = load_data([Analysis_Path sigma4{i} '-cf/' sigma4{i} '-cf']);
            fig_offset = 3;
        end
        
        if i == 1
            subplot_offset = 0;
            cells = [1922 3181 5793 7366];
            type = [1 2 1 2];
        else
            subplot_offset = 4;
            cells = [320 2461 3377 4411];% 5641];
            type = [1 2 1 2 2];
        end
        
        datarun = load_neurons(datarun);
        datarun_cf = load_neurons(datarun_cf);
        sigma_data = interleaved_data_prep(datarun, 1100, 120, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);
        sigma_data_cf = interleaved_data_prep(datarun_cf, 1100, 120, 'cell_spec', 'all', 'visual_check', 0);
        
        for i_cell = 1:length(cells)
            try
                % replace with cell findered data
                sigma_data.testspikes(:,i_cell) = sigma_data_cf.testspikes(:,cells(i_cell) == datarun_cf.cell_ids);
            catch
                disp(['cell ' num2str(cells(i_cell)) ' not cell findered'])
            end
        end
        clear sigma_data_cf
        
        % separate the conditions
        mask_data.testspikes = sigma_data.testspikes(1:30, :);
        LES_data.testspikes = sigma_data.testspikes(31:60, :);
        comp_data.testspikes = sigma_data.testspikes(61:90, :);
        compLES_data.testspikes = sigma_data.testspikes(91:end, :);
        
        % reg data
        %         datarun = load_data([ Analysis_Path regA{i} '/' regA{i}]);
        %         datarun = load_neurons(datarun);
        %         regA_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells{i}, 'datarun_class', class_datarun, 'visual_check', 0);
        datarun = load_data([ Analysis_Path regB{i} '/' regB{i}]);
        datarun_cf = load_data([ Analysis_Path regB{i} '-cf/' regB{i} '-cf']);
        datarun = load_neurons(datarun);
        datarun_cf = load_neurons(datarun_cf);
        regB_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);
        regB_data_cf = interleaved_data_prep(datarun_cf, 1100, 30, 'cell_spec', 'all', 'visual_check', 0);
        for i_cell = 1:length(cells)
            try
                % replace with cell findered data
                regB_data.testspikes(:,i_cell) = regB_data_cf.testspikes(:,cells(i_cell) == datarun_cf.cell_ids);
            catch
                disp(['cell ' num2str(cells(i_cell)) ' not cell findered sigma'])
            end
        end
        clear regB_data_cf
        
        % plot for each cell
        for cell = 1%:size(sigma_data.testspikes,2)
            
            % reg versus mask
            figure(1);
            reg = IDP_plot_PSTH(regB_data, cell);
            mask = IDP_plot_PSTH(mask_data, cell);
            comp = IDP_plot_PSTH(comp_data, cell);
            LES = IDP_plot_PSTH(LES_data, cell);
            compLES = IDP_plot_PSTH(compLES_data, cell);
            
            figure(2); hold on; plot(IDP_spike_counts(reg, image_transitions),IDP_spike_counts(mask, image_transitions), '.', 'Color', default_colors(type(cell), :))
            figure(3); hold on; plot(IDP_spike_counts(reg, image_transitions),IDP_spike_counts(mask+comp, image_transitions), '.', 'Color', default_colors(type(cell), :))
            
            
            if 0
                
                % reg versus mask
                figure(1);
                plot(reg)
                hold on
                plot(mask)
                legend('reg', 'mask')
                exportfig(gcf, [fig_save '/maskin' num2str(mod(i,2)+1) '_sigma' num2str(2*sigma) '_regmask_cell_' num2str(cell)], 'Bounds', 'loose', 'Color', 'rgb')
                
                
                % reg versus complement + mask
                figure(2);
                plot(reg,'Color', default_colors(4,:))
                hold on
                plot(mask+comp, 'Color', default_colors(5,:))
                plot(-mask, 'Color', default_colors(5,:)+[0.0 0.0 0.3])
                plot(-comp, 'Color', default_colors(5,:)-[0.3 0.3 0.15])
                hold off
                legend('Reg', 'Mask+Comp', 'Mask', 'Comp')
                exportfig(gcf, [fig_save '/maskin' num2str(mod(i,2)+1) '_sigma' num2str(2*sigma) '_compmask_cell_' num2str(cell)], 'Bounds', 'loose', 'Color', 'rgb')
                
                % mask versus LES
                figure(3); hold on;
                plot(mask);
                plot(LES);
                legend('mask', 'LES')
                exportfig(gcf, [fig_save '/maskin' num2str(mod(i,2)+1) '_sigma' num2str(2*sigma) '_maskLES_cell_' num2str(cell)], 'Bounds', 'loose', 'Color', 'rgb')
                
                
                % reg versus LES+complement
                figure(4); hold on;
                plot(reg);
                plot(compLES);
                legend('reg', 'compLES')
                exportfig(gcf, [fig_save '/maskin' num2str(mod(i,2)+1) '_sigma' num2str(2*sigma) '_regcompLES_cell_' num2str(cell)], 'Bounds', 'loose', 'Color', 'rgb')
                
                
                figure(5); hold on;
                plot(compLES,'Color', default_colors(4,:));
                plot(comp+LES, 'Color', default_colors(5,:));
                plot(-LES, 'Color', default_colors(5,:)+[0.0 0.0 0.3])
                plot(-comp, 'Color', default_colors(5,:)-[0.3 0.3 0.15])
                legend('compLES', 'LES+Comp', 'LES', 'Comp')
                exportfig(gcf, [fig_save '/maskin' num2str(mod(i,2)+1) '_sigma' num2str(2*sigma) '_compLES_cell_' num2str(cell)], 'Bounds', 'loose', 'Color', 'rgb')
                
                close all
            end
            
        end
        figure(2)
        plot([0 2000], [0 2000])
        axis([0 2000 0 2000])
        axis square
        set(gcf, 'Position', [2 2 2 2])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
        figure(3)
        plot([0 2000], [0 2000])
        axis([0 2000 0 2000])
        axis square
        set(gcf, 'Position', [2 2 2 2])
        set(gca, 'XTick', [])
        set(gca, 'YTick', [])
    end
end

%% look at parts of the movie
% first load the movie and the image transitions
default_colors = get(gca,'ColorOrder');
type_str = {'ON', 'OFF'};
for sigma = 1%:2
    % load masks
    if sigma == 1
        if maskin == 1
            mask_file = load('/Users/Nora/Desktop/Data/2016-01-05-0-Visual/2016-01-05-0/testmask_1922_stix2/mask_and_filter.mat');
        else
            mask_file = load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/mask_320_sigma2/mask_and_filter.mat');
        end
    else
        if maskin == 1
            mask_file = load('/Users/Nora/Desktop/Data/2016-01-05-0-Visual/2016-01-05-0/testmask_1922_sigma4/mask_and_filter.mat');
        else
            mask_file = load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/mask_320_sigma4/mask_and_filter.mat');
        end
    end
    % Load the data
    for i = 1%:2
        if sigma==1
            datarun = load_data([Analysis_Path sigma2{i} '/' sigma2{i}]);
            datarun_cf = load_data([Analysis_Path sigma2{i} '-cf/' sigma2{i} '-cf']);
            fig_offset = 0;
        else
            datarun = load_data([Analysis_Path sigma4{i} '/' sigma4{i}]);
            datarun_cf = load_data([Analysis_Path sigma4{i} '-cf/' sigma4{i} '-cf']);
            fig_offset = 3;
        end
        if i == 1
            subplot_offset = 0;
            cells = [1922 3181 5793 7366];
            type = [1 2 1 2];
        else
            subplot_offset = 4;
            cells = [320 2461 3377 4411];% 5641];
            type = [1 2 1 2 2];
        end
        datarun = load_neurons(datarun);
        datarun_cf = load_neurons(datarun_cf);
        sigma_data = interleaved_data_prep(datarun, 1100, 120, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);
        sigma_data_cf = interleaved_data_prep(datarun_cf, 1100, 120, 'cell_spec', 'all', 'visual_check', 0);
        for i_cell = 1:length(cells)
            try
                % replace with cell findered data
                sigma_data.testspikes(:,i_cell) = sigma_data_cf.testspikes(:,cells(i_cell) == datarun_cf.cell_ids);
            catch
                disp(['cell ' num2str(cells(i_cell)) ' not cell findered'])
            end
        end
        clear sigma_data_cf
        % separate the conditions
        mask_data.testspikes = sigma_data.testspikes(1:30, :);
        LES_data.testspikes = sigma_data.testspikes(31:60, :);
        comp_data.testspikes = sigma_data.testspikes(61:90, :);
        compLES_data.testspikes = sigma_data.testspikes(91:end, :);
        % reg data
        datarun = load_data([ Analysis_Path regB{i} '/' regB{i}]);
        datarun_cf = load_data([ Analysis_Path regB{i} '-cf/' regB{i} '-cf']);
        datarun = load_neurons(datarun);
        datarun_cf = load_neurons(datarun_cf);
        regB_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);
        regB_data_cf = interleaved_data_prep(datarun_cf, 1100, 30, 'cell_spec', 'all', 'visual_check', 0);
        for i_cell = 1:length(cells)
            try
                % replace with cell findered data
                regB_data.testspikes(:,i_cell) = regB_data_cf.testspikes(:,cells(i_cell) == datarun_cf.cell_ids);
            catch
                disp(['cell ' num2str(cells(i_cell)) ' not cell findered sigma'])
            end
        end
        clear regB_data_cf
        
        % plot for each cell
        for cell = 1:length(cells)
            
            % reg versus mask
            figure(1);
            reg = IDP_plot_PSTH(regB_data, cell);
            hold on
            mask = IDP_plot_PSTH(mask_data, cell);
            comp = IDP_plot_PSTH(comp_data, cell);
            for i_line = image_transitions
                plot([i_line i_line], [0 200], 'Color', default_colors(4,:))
            end
            legend('reg', 'mask', 'comp')
            [times, ~] = ginput();
            
            for i_time = 1:length(times)
                figure;
                current_image_begin = image_transitions(find(diff(ceil(image_transitions/times(i_time)))==1, 1)); 
                subplot(2,1,1);imagesc(NSmovie(:,:,current_image_begin-1)+500*edge(mask_file.mask(:,:,cell)));caxis([0 255]); axis image;colormap gray; title('Previous Image');
                subplot(2,1,2);imagesc(NSmovie(:,:,current_image_begin)+500*edge(mask_file.mask(:,:,cell)));caxis([0 255]); axis image;colormap gray; title('Current Image');
                time_from_saccade = num2str(times(i_time) - current_image_begin);
                bin = round(times(i_time));
                title_str = [type_str{type(cell)} '; Frames from Saccade ' time_from_saccade '; Reg' num2str(reg(bin)), 'Comp', num2str(comp(bin)), 'Mask', num2str(mask(bin))];
                suptitle(title_str);
                set(gcf, 'Position', [2 2 8 16])
            end
            pause()
            close all
            
        end
    end
end



%% look at amacrine cells

% first load the movie and the image transitions
default_colors = get(gca,'ColorOrder');
type_str = {'ON', 'OFF'};
for sigma = 1%:2
    % load masks
    if sigma == 1
        if maskin == 1
            mask_file = load('/Users/Nora/Desktop/Data/2016-01-05-0-Visual/2016-01-05-0/testmask_1922_stix2/mask_and_filter.mat');
        else
            mask_file = load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/mask_320_sigma2/mask_and_filter.mat');
        end
    else
        if maskin == 1
            mask_file = load('/Users/Nora/Desktop/Data/2016-01-05-0-Visual/2016-01-05-0/testmask_1922_sigma4/mask_and_filter.mat');
        else
            mask_file = load('/Volumes/Data/2016-01-05-0/Visual/2016-01-05-0/mask_320_sigma4/mask_and_filter.mat');
        end
    end
    % Load the data
    for i = 1%:2
        if sigma==1
            datarun = load_data([Analysis_Path sigma2{i} '/' sigma2{i}]);
            datarun_cf = load_data([Analysis_Path sigma2{i} '-cf/' sigma2{i} '-cf']);
            fig_offset = 0;
        else
            datarun = load_data([Analysis_Path sigma4{i} '/' sigma4{i}]);
            datarun_cf = load_data([Analysis_Path sigma4{i} '-cf/' sigma4{i} '-cf']);
            fig_offset = 3;
        end
        if i == 1
            subplot_offset = 0;
            cells = [1922 3181 5793 7366];
            type = [1 2 1 2];
        else
            subplot_offset = 4;
            cells = [320 2461 3377 4411];% 5641];
            type = [1 2 1 2 2];
        end
        datarun = load_neurons(datarun);
        datarun_cf = load_neurons(datarun_cf);
        sigma_data = interleaved_data_prep(datarun, 1100, 120, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);
        sigma_data_cf = interleaved_data_prep(datarun_cf, 1100, 120, 'cell_spec', 'all', 'visual_check', 0);
        for i_cell = 1:length(cells)
            try
                % replace with cell findered data
                sigma_data.testspikes(:,i_cell) = sigma_data_cf.testspikes(:,cells(i_cell) == datarun_cf.cell_ids);
            catch
                disp(['cell ' num2str(cells(i_cell)) ' not cell findered'])
            end
        end
        clear sigma_data_cf
        % separate the conditions
        mask_data.testspikes = sigma_data.testspikes(1:30, :);
        LES_data.testspikes = sigma_data.testspikes(31:60, :);
        comp_data.testspikes = sigma_data.testspikes(61:90, :);
        compLES_data.testspikes = sigma_data.testspikes(91:end, :);
        % reg data
        datarun = load_data([ Analysis_Path regB{i} '/' regB{i}]);
        datarun_cf = load_data([ Analysis_Path regB{i} '-cf/' regB{i} '-cf']);
        datarun = load_neurons(datarun);
        datarun_cf = load_neurons(datarun_cf);
        regB_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells, 'datarun_class', class_datarun, 'visual_check', 0);
        regB_data_cf = interleaved_data_prep(datarun_cf, 1100, 30, 'cell_spec', 'all', 'visual_check', 0);
        for i_cell = 1:length(cells)
            try
                % replace with cell findered data
                regB_data.testspikes(:,i_cell) = regB_data_cf.testspikes(:,cells(i_cell) == datarun_cf.cell_ids);
            catch
                disp(['cell ' num2str(cells(i_cell)) ' not cell findered sigma'])
            end
        end
        amacrine_reg_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', 'Off Amacrine', 'datarun_class', class_datarun, 'visual_check', 0);
        amacrine_avg.testspikes{1} = cell2mat(amacrine_reg_data.testspikes(:));
        amacrine = IDP_plot_PSTH(amacrine_avg, 1);
        clear regB_data_cf
        
        % plot for each cell
        for cell = 1:length(cells)
            
            % reg versus mask
            figure(1);
            reg = IDP_plot_PSTH(regB_data, cell);
            hold on;
            mask = IDP_plot_PSTH(mask_data, cell); hold on;
            comp = IDP_plot_PSTH(comp_data, cell); hold on;
            plot(amacrine/400);
            for i_line = image_transitions
                plot([i_line i_line], [0 200], 'Color', default_colors(5,:))
            end
            legend('reg', 'mask', 'comp')
            %[times, ~] = ginput();
            pause()
            close all
            
        end
    end
end

 
%%
%% look at all the rasters!

% first 30 are mask
% next 30 are LES
% next 30 are complement
% final 30 are complement + LES
% Load the data
for i = 1:2
    for sigma = 1:2
        if sigma==1
            datarun = load_data([Analysis_Path sigma2{i} '/' sigma2{i}]);
            fig_offset = 0;
        else
            datarun = load_data([Analysis_Path sigma4{i} '/' sigma4{i}]);
            fig_offset = 3;
        end
        
        datarun = load_neurons(datarun);
        sigma_data = interleaved_data_prep(datarun, 1100, 120, 'cell_spec', cells{i}, 'datarun_class', class_datarun, 'visual_check', 0);
        clear datarun
    % separate the conditions
    mask_data{sigma}.testspikes = sigma_data.testspikes(1:30, :);
    %LES_data.testspikes = sigma_data.testspikes(31:60, :);
    comp_data{sigma}.testspikes = sigma_data.testspikes(61:90, :);
    %compLES_data.testspikes = sigma_data.testspikes(91:end, :);
    end
    
    % reg data
    %         datarun = load_data([ Analysis_Path regA{i} '/' regA{i}]);
    %         datarun = load_neurons(datarun);
    %         regA_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells{i}, 'datarun_class', class_datarun, 'visual_check', 0);
    datarun = load_data([ Analysis_Path regB{i} '/' regB{i}]);
    datarun = load_neurons(datarun);
    regB_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cells{i}, 'datarun_class', class_datarun, 'visual_check', 0);
    
    % plot for each cell
    for cell = 1:size(mask_data{1}.testspikes,2)
        figure(2); hold on 
        figure(1); hold on
        reg = IDP_plot_PSTH(regB_data, cell, 'color', 0);
        for sigma=1:2
           color = [1 1 1]*(1-sigma/12)-0.2;
           figure(1); mask = IDP_plot_PSTH(mask_data{sigma}, cell,'color',color);
           NSEM_Corr(cell, sigma) = err(reg, mask);
           figure(2); comp = IDP_plot_PSTH(comp_data{sigma}, cell,'color', color);
           surround_struct(cell, sigma) = std(comp);
        end
        figure(1);
        reg = IDP_plot_PSTH(regB_data, cell, 'color', [0 0 0]);
        exportfig(gcf, [fig_save '/maskin' num2str(mod(i,2)+1) '_mask_cell_' num2str(cell)], 'Bounds', 'loose', 'Color', 'rgb');
        clf;
        figure(2);
        exportfig(gcf, [fig_save '/maskin' num2str(mod(i,2)+1) '_comp_cell_' num2str(cell)], 'Bounds', 'loose', 'Color', 'rgb');
        clf;
        
    end
    save(['SpotStats2016-01-05-0' num2str(i)], 'NSEM_Corr', 'surround_struct');
end




%%

for fig = 1%:10
    figure(fig);
    % exportfig(gcf, [fig_save '/fig' num2str(fig)], 'Bounds', 'loose', 'Color', 'rgb')
end

