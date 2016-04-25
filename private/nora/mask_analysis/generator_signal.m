clear 
set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 12 4])
GS = 1;
sta_scale = 4;
colors = colormap(winter(6));

TCF = [];
TMG = [];
TMF = [];
TCG = [];

params_list

for exp = exps
    disp(exp{1})
    eval(exp{1})

n_reg = length(reg);
mkdir(fig_save);
n_subgroups = length(cell_idx);



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


%% load movie
if GS
    i_chunk = 1;
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
    NSmovie = movie_chunk;
    interval(i_chunk) = size(movie_chunk,3);
    i_chunk = 2;
    
    % mask movie
    while size(NSmovie,3) < 1200
        load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
        interval(i_chunk) = size(movie_chunk,3);
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

interval_frame = cumsum(interval);

%%
% plot the PSTH for each condition for each cell

for subgroup = 1:n_subgroups
    for i_cell = cell_idx{subgroup}
        disp(cells_orig(i_cell))
        master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
        [reg, time] = IDP_plot_PSTH(reg_data{1},i_cell, 'color', 0, 'smoothing', 10);
        %[params,~,~] = fit_sta(datarun_class.stas.stas{master_idx}, 'fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', true);
        %sta = squeeze(sum(datarun_class.stas.stas{master_idx},3));
        %%{
        sta_fit = datarun_class.vision.sta_fits{master_idx}; 
        ideal_center = make_Gaussian_two_d('center_point_x', sta_fit.mean(1), 'center_point_y', 40 -sta_fit.mean(2), 'sd_x', sta_fit.sd(1), 'sd_y', sta_fit.sd(2), 'rotation_angle', sta_fit.angle, 'x_dim', 80, 'y_dim', 40);
        ideal_surround = make_Gaussian_two_d('center_point_x', sta_fit.mean(1), 'center_point_y', 40-sta_fit.mean(2), 'sd_x', 2*sta_fit.sd(1), 'sd_y',2*sta_fit.sd(2), 'rotation_angle', sta_fit.angle, 'x_dim', 80, 'y_dim', 40);
        time_course=mean([datarun_class.vision.timecourses(master_idx).r ...
            datarun_class.vision.timecourses(master_idx).g ...
            datarun_class.vision.timecourses(master_idx).b],2);
        sta_center = reshape((ideal_center(:))*time_course', [40 80 30]);
        sta_surround = reshape((ideal_surround(:))*[time_course(3:end);0; 0]', [40 80 30]);
        %}
        sta_center = flip(sta_center,1);  
        sta_center = flip(sta_center,2); 
        sta_center = flip(sta_center,3); 
        sta_surround = flip(sta_surround,1);  
        sta_surround = flip(sta_surround,2); 
        sta_surround = flip(sta_surround,3);
        for i = 1%:length(mask_conditions{subgroup})
            s = sigmas(mask_conditions{subgroup}(i));
            if strcmp(piece, '2016-04-21-1')
                if s<8
                    load(['/Volumes/Data/2016-04-21-1/Visual/2016-04-21-1_NJB_Masks/Maskin_allcells_sigma' num2str(s) '.mat'])
                else
                    load(['/Volumes/Data/2016-04-21-1/Visual/2016-04-21-1_NJB_Masks/Maskin_cells' num2str(subgroup) '_sigma' num2str(s) '.mat'])
                end
            else
                if s==2
                    load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_allcells_sigma' num2str(s) '.mat'])
                else
                    load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_cells' num2str(subgroup) '_sigma' num2str(s) '.mat']);
                end
            end
            
            comp_mask = mod(mask+1,2);
            mask = imresize(mask, 1/sta_scale, 'box');
            mask = repmat(mask, 1, 1, 1200);
            comp_mask = imresize(comp_mask, 1/sta_scale, 'box');
            comp_mask = repmat(comp_mask, 1, 1, 1200);
            mask_movie = NSmovie .* mask;
            mask_gen_signal_center = squeeze(convn(mask_movie, sta_center, 'valid'));
            mask_gen_signal_surround = squeeze(convn(mask_movie, sta_surround, 'valid'));
            reg_gen_signal_center = squeeze(convn(NSmovie, sta_center, 'valid')); 
            reg_gen_signal_surround = squeeze(convn(NSmovie, sta_surround, 'valid')); 
            comp_gen_signal_center = squeeze(convn(NSmovie.*comp_mask, sta_center, 'valid')); 
            comp_gen_signal_surround = squeeze(convn(NSmovie.*comp_mask, sta_surround, 'valid')); 
            mask_PSTH = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
            comp = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
            
            close all
                        
            model_arch = 'y~log(1+exp(b1*x1+b2*x2))';
            init = [6 -1.5]/4;
            responses = reg(30:end);
            data_matrix = [reg_gen_signal_center(1:1166) reg_gen_signal_surround(1:1166)];
            model = fitnlm(data_matrix, responses, model_arch, init);
            %NLN_model = fitnlm(data_matrix, responses, 'y~log(1+exp(b1*x1))+log(1+exp(b2*x2))', [6 -1.5]);
            %figure; hold on; plot([responses, predict(model, data_matrix)]);
            %Coeff(i_cell, i,1, :) = model.Coefficients.Estimate;
            Corr(i_cell, i,1, :) = model.Rsquared.Ordinary;
            reg_gen_signal = data_matrix*model.Coefficients.Estimate;
            
            responses = squeeze(mask_PSTH(30:end));
            data_matrix = [mask_gen_signal_center(1:1166) mask_gen_signal_surround(1:1166)];
            model = fitnlm(data_matrix, responses, model_arch, init);
            %figure; plot([responses, predict(model, data_matrix)]);
            %Coeff(i_cell, i,2, :) = model.Coefficients.Estimate;
            Corr(i_cell, i,2, :) = model.Rsquared.Ordinary;
            mask_gen_signal = data_matrix*model.Coefficients.Estimate;
            
            responses = squeeze(comp(30:end));
            data_matrix = [comp_gen_signal_center(1:1166) comp_gen_signal_surround(1:1166)];
            model = fitnlm(data_matrix, responses, model_arch, init);
            %figure; plot([responses, predict(model, data_matrix)]);
            %Coeff(i_cell, i,3, :) = model.Coefficients.Estimate;
            Corr(i_cell, i,3, :) = model.Rsquared.Ordinary;
            comp_gen_signal = data_matrix*model.Coefficients.Estimate;
            close all
           
            
            start = 1;
            for i_interval = 1:24
               comp_firing(i_interval) = sum(comp(start:interval_frame(i_interval)-1));
               mask_firing(i_interval) = sum(mask(start:interval_frame(i_interval)-1));
               try
                   mask_gen(i_interval) = sum(mask_gen_signal((start:interval_frame(i_interval)-1)-30));
                   comp_gen(i_interval) = sum(comp_gen_signal((start:interval_frame(i_interval)-1)-30));
               catch 
                   mask_gen(i_interval) = 0;
                   comp_gen(i_interval) = 0;
               end
               start = interval_frame(i_interval);
            end
            TCF = [TCF; comp_firing];
            TMF = [TMF; mask_firing];
            TMG = [TMG; mask_gen];
            TCG = [TCG; comp_gen];
            
            %{
            plot(time(30:end),comp(30:end))
            hold on;
            plot(time(30:end), mask_gen_signal(1:1166))
            xlim([1 9])
            axis off
            exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) 'sigma' num2str(i) '_compVSmaskGS'], 'Bounds', 'loose', 'Color', 'rgb')
            close all
            
            plot(time(30:end),reg(30:end))
            hold on;
            plot(time(30:end),mask_PSTH(30:end))
            plot(time(30:end), comp_gen_signal(1:1166))
            xlim([1 9])
            axis off
            exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) 'sigma' num2str(i) '_maskVScompGS'], 'Bounds', 'loose', 'Color', 'rgb')
            %}
            
        end
        
    end
end
clear reg 
end

plot(TCF, TMG, '.k')
