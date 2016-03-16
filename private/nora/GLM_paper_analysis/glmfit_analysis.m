dsave = '/Volumes/Lab/Users/Nora/GLMFits/';
%{
piece = '2016-02-17-1/';
Analysis_Path = '/Volumes/Analysis/2016-02-17-1/data022-data028';
datarun_class = load_data([Analysis_Path '/data024/data024'], struct('load_neurons', 0, 'load_params', 1));
%}
%%{
piece = '2015-05-27-3/';
Analysis_Path = '/Volumes/Analysis/2015-05-27-3/data001-data005';
datarun_class = load_data([Analysis_Path '/data001/data001'], struct('load_neurons', 0, 'load_params', 1));
 %}

for cell_type = {'On Parasol'}
    cells = get_cell_ids(datarun_class, cell_type);
    n_cells = length(cells);
    BPS = zeros(2, n_cells);
    for i_cell = 1:n_cells
       load([dsave piece num2str(cells(i_cell)) '.mat']);
       %[~, ~, crm_bps, ~] = rastbps_comp_findPS(fittedGLM.xvalperformance.rasters.recorded,fittedGLM.t_bin,fittedGLM.rawfit.ps_basis); 
       BPS(1,i_cell) = fittedGLM.xvalperformance.logprob_glm_bpspike;%/crm_bps;
       load([dsave piece num2str(cells(i_cell)) 'NSEM.mat']);
       %[~, ~, crm_bps, ~] = rastbps_comp_findPS(fittedGLM.xvalperformance.rasters.recorded,fittedGLM.t_bin,fittedGLM.rawfit.ps_basis); 
       BPS(2,i_cell) = fittedGLM.xvalperformance.logprob_glm_bpspike;%/crm_bps;
    end
    figure; plot(BPS(1,:), BPS(2,:), '.')
    save([dsave piece 'BPSOn.mat'], 'BPS')
end


