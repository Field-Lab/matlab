clear
set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 6 4])
sta_scale = 4;

%%
params_list
for exp = [2 3 5]
    eval(exps{exp}) % 2 3 5
    n_reg = length(reg);
    mkdir(fig_save);
    n_subgroups = length(cell_idx);
    
    %%
    % load classification run
    datarun_class = load_data(class);
    datarun_class = load_params(datarun_class);
    datarun_class = load_sta(datarun_class, 'load_sta', cells_orig, 'save_rf', 1);

    %%
    % plot the PSTH for each condition for each cell
    
    glm_files = dir(['/Volumes/Lab/Users/Nora/GLMFits_masking/' piece '/NSEM_full_opt/*NSEM_optmodel.mat']);
    if strcmp(piece, '2016-04-21-1')
        glm_fit_idx = [6 3 5 1 2 4];
    end
    %%
    % 1 = center, 2 = surround
    % 1 = reg, 2 = spot, 3 = gap
    for subgroup = 1:n_subgroups
        for i_cell = cell_idx{subgroup}
            disp(cells_orig(i_cell))
               master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
               if strcmp(piece, '2016-04-21-1')
                    load(['/Volumes/Lab/Users/Nora/GLMFits_masking/2016-04-21-1/NSEM_full_opt/' glm_files(glm_fit_idx(i_cell)).name]);
                else
                    load(['/Volumes/Lab/Users/Nora/GLMFits_masking/' piece '/NSEM_full_opt/' num2str(cells_orig(i_cell)) 'NSEM_optmodel.mat']);
                end
                sta_fit = opt_model.sta_fit.params;
                plot(1000*(0:(1/119.5):(29/119.5)),opt_model.new_time,'k')
                set(gca, 'YTick', []);
                %exportfig(gcf, '/Volumes/Lab/Users/Nora/time_filter.eps', 'Bounds', 'loose', 'Color', 'rgb')
                imagesc(opt_model.orig_glm.linearfilters.Stimulus.space_rk1)
                colormap parula
                axis([10 30 2 22])
                axis square
                axis off
                %exportfig(gcf, '/Volumes/Lab/Users/Nora/space_filter.eps', 'Bounds', 'tight', 'Color', 'rgb', 'Renderer', 'zbuffer')

        end
    end
end