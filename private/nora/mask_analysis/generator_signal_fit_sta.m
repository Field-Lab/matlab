clear 
params_201602176
n_reg = length(reg);
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

n_subgroups = length(cell_idx);
GS = 1;
sta_scale = 4;
default_colors = get(gca,'ColorOrder');


%% load and organize repeats 

% load reg repeats
for i_reg = 1:n_reg
    datarun{i_reg} = load_data([ Analysis_Path reg{i_reg} '/' reg{i_reg}]);
    datarun{i_reg} = load_neurons(datarun{i_reg});
    reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells_reg{i_reg}, 'visual_check', 1);
    clear datarun
end

% load masking repeats
datarun = load_data([ Analysis_Path masking '/' masking]);
datarun = load_neurons(datarun);
mask_data = interleaved_data_prep(datarun, 1100, n_masks*2*30, 'cell_spec', cells_masking, 'visual_check', 0);
clear datarun
% separate each condition
idx = 1:30;
for count = 1:(n_masks*2)
    condition{count}.testspikes = mask_data.testspikes(idx,:); idx = idx+30;
end

%%
% load classification run 
datarun_class = load_data(class);
datarun_class = load_params(datarun_class);
% avg RF
datarun_class = load_sta(datarun_class);

%%
%{
avg_rf = sum(get_average_rf(datarun_class, 'On Parasol', 'scale', 5),3);

n_angles = 10;
slice = zeros(201,1);
range = -100:100;

for i = 1:n_angles
    figure(1); imagesc(avg_rf); axis image; title(i)
    prof = improfile;
    prof = [zeros(100,1); prof; zeros(100,1)];
    [~,center] = max(prof);
    slice = slice + prof(center+range);
end
slice = slice/n_angles;
plot(slice)

a = fit((1:201)'/5,slice, fittype('gauss2'));
width = a.c1/sqrt(2);
avg_profile = [slice, (1:201)'/(5*width)];
%}

%% load movie
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
    NSmovie = NSmovie - 64;
    NSmovie = imresize(NSmovie, 1/sta_scale, 'box');
    NSmovie = permute(NSmovie, [2 1 3]);
    
end


%%
for subgroup = 1:n_subgroups
    for i_cell = cell_idx{subgroup}
        disp(cells_orig(i_cell))
        master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
        [reg, time] = IDP_plot_PSTH(reg_data{1},i_cell, 'color', 0, 'smoothing', 10);
        [params,sta,~] = fit_sta(datarun_class.stas.stas{master_idx}, 'fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', true);
        sta = squeeze(sum(sta,3));
        sta = flip(sta,1);  
        sta = flip(sta,2); 
        sta = flip(sta,3);
        for i = 1:length(mask_conditions{subgroup})
            s = sigmas(mask_conditions{subgroup}(i));
            if s==2
                load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_allcells_sigma' num2str(s) '.mat'])
            else
                load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_cells' num2str(subgroup) '_sigma' num2str(s) '.mat']);
            end
            comp_mask = mod(mask+1,2);  
            mask = imresize(mask, 1/sta_scale, 'box');
            mask = repmat(mask, 1, 1, 1200);
            comp_mask = imresize(comp_mask, 1/sta_scale, 'box');
            comp_mask = repmat(comp_mask, 1, 1, 1200);
            mask_movie = NSmovie .* mask;
            mask_gen_signal_center = squeeze(convn(mask_movie, sta, 'valid')); 
            reg_gen_signal_center = squeeze(convn(NSmovie, sta, 'valid')); 
            comp_gen_signal_center = squeeze(convn(NSmovie.*comp_mask, sta, 'valid')); 
            
            %comp_gen_signal_center = mean(reshape(NSmovie.*comp_mask, 40*80, 1200),1); 
            mask_PSTH = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
            comp = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
            
            model = fitnlm(reg_gen_signal_center(1:1166), reg(30:end), 'y~b1*log(1+exp(b2*x1-b3))', [1 1 1]);

        end
    end
end 