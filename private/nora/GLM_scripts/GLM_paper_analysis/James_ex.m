%%
fit_dir = '/Volumes/Lab/Users/akheitman/NSEM_Home/GLM_Output_Analysis/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/';
folders = dir([fit_dir 'NSEM_mapPRJ/']);
default_colors = get(gca,'ColorOrder');
hold on;
for i_exp = 1:length(folders)
    exp = folders(i_exp).name;
    if ~strcmp(exp(1), '.')
        for type = {'ON', 'OFF'}
            cells = dir([fit_dir 'NSEM_mapPRJ/' exp '/' type{1} '*.mat']);
            n_cells = length(cells);
            BPS = zeros(2, n_cells);
            for i_cell = 1:n_cells
                disp(i_cell)
                load([fit_dir 'WN_mapPRJ/' exp '/' cells(i_cell).name]);
                [~, ~, crm_bps, ~] = rastbps_comp_findPS(fittedGLM.xvalperformance.rasters.recorded,fittedGLM.t_bin,fittedGLM.rawfit.ps_basis);
                BPS(1,i_cell) = fittedGLM.xvalperformance.logprob_glm_bpspike/crm_bps;
                load([fit_dir 'NSEM_mapPRJ/' exp '/' cells(i_cell).name]);
                [~, ~, crm_bps, ~] = rastbps_comp_findPS(fittedGLM.xvalperformance.rasters.recorded,fittedGLM.t_bin,fittedGLM.rawfit.ps_basis);
                BPS(2,i_cell) = fittedGLM.xvalperformance.logprob_glm_bpspike/crm_bps;
            end
            if strcmp(type, 'ON');
                plot(BPS(1,:), BPS(2,:), '.', 'Color', default_colors(i_exp,:));
            else
                plot(BPS(1,:), BPS(2,:), '+', 'Color', default_colors(i_exp,:));
            end
        end
    end
end
    
%%
dsave = '/Volumes/Lab/Users/Nora/GLMFits/';
%%{
piece = '2016-02-17-1/';
Analysis_Path = '/Volumes/Analysis/2016-02-17-1/data022-data028';
datarun_class = load_data([Analysis_Path '/data024/data024'], struct('load_neurons', 0, 'load_params', 1));
%}

for cell_type = {'Off Parasol'}
    cells = get_cell_ids(datarun_class, cell_type);
    n_cells = length(cells);
    BPS = zeros(2, n_cells);
    for i_cell = 2:n_cells
        disp(i_cell)
       load([dsave piece num2str(cells(i_cell)) '.mat']);
try
       [~, ~, crm_bps, ~] = rastbps_comp_findPS(fittedGLM.xvalperformance.rasters.recorded,fittedGLM.t_bin,fittedGLM.rawfit.ps_basis); 
       BPS(1,i_cell) = fittedGLM.xvalperformance.logprob_glm_bpspike;%/crm_bps;
%       load([dsave piece num2str(cells(i_cell)) 'NSEM.mat']);
%       [~, ~, crm_bps, ~] = rastbps_comp_findPS(fittedGLM.xvalperformance.rasters.recorded,fittedGLM.t_bin,fittedGLM.rawfit.ps_basis); 
       BPS(2,i_cell) = fittedGLM.xvalperformance.logprob_glm_bpspike;%/crm_bps;
catch
     disp([num2str(i_cell) ' crashed'])   
end
end
    figure; plot(BPS(1,:), BPS(2,:), '.')
    %save([dsave piece 'BPSOff.mat'], 'BPS')
end


