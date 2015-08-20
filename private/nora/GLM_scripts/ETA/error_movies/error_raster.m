function [testmovie, res] = error_raster(exp,celltype)

stimtype = 'NSEM';

datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/NSEM_mapPRJ/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];

% Get all fits of that cell type
matfiles=dir([datapath exp_names(exp,:) celltype '*.mat']);
n_cells=length(matfiles);

% Load one fittedGLM to get the experiment info
load([datapath exp_names(exp,:) matfiles(1).name]);
GLMType = fittedGLM.GLMType;
exp_nm = fittedGLM.cellinfo.exp_nm;

% Load and process stimulus
[StimulusPars, ~] = StimulusParams(exp_nm, stimtype, GLMType.map_type);
[movie, ~, ~]          = loadmoviematfile(exp_nm , stimtype, GLMType.cone_model,'testmovie');
testmovie             = movie{1}.matrix(:,:,StimulusPars.slv.testframes);

% Loop through the cells and find the residual spikes for each cell.
res.spikes = zeros(n_cells, 10*length(testmovie));
res.cells = zeros(n_cells, 1);
res.centers = zeros(n_cells, 2);
for i_cell = 1:n_cells
    disp(i_cell)
    % Load fittedGLM
    load([datapath exp_names(exp,:) matfiles(i_cell).name]);
    % Calculate the rate from GLM
    res.spikes(i_cell,:) = glm_error(fittedGLM);
    res.cells(i_cell) = fittedGLM.cellinfo.cid;
    res.centers(i_cell,1) = fittedGLM.cellinfo.slave_centercoord.y_coord;
    res.centers(i_cell,2) = fittedGLM.cellinfo.slave_centercoord.x_coord;
end

end
