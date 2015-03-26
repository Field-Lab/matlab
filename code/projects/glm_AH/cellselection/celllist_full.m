% Junk Script


%%% Purpose  %%%
% Load an easy to access allcells data structure
% So we don't need to call load_data to access cell ids f

% AKHEITMAN 2014-11-11



BD = NSEM_BaseDirectories;

savedir = sprintf('%s',BD.Cell_Selection)


%
allcells = cell(4,1);
cellselectiontype = 'all';
for i_exp = 1:4
    
    [exp_nm,~,expname,badcells]  = cell_list( i_exp, cellselectiontype);  
    allcells{i_exp,1}.exp_nm = exp_nm;
    allcells{i_exp,1}.expname = expname;
    
    [~, DirPars_WN   datarun_slv_WN   datarun_mas] = Directories_Params_v23_Conv(exp_nm, 'WN', 'mapPRJ',1);
    [~, DirPars_NSEM datarun_slv_NSEM datarun_mas] = Directories_Params_v23_Conv(exp_nm, 'NSEM', 'mapPRJ',1);
    ONP  = intersect(datarun_slv_WN.cell_types{1}.cell_ids , datarun_slv_NSEM.cell_types{1}.cell_ids);
    OFFP = intersect(datarun_slv_WN.cell_types{2}.cell_ids , datarun_slv_NSEM.cell_types{2}.cell_ids);
    for i_cell = 1:length(badcells);
        ONP(find(ONP == badcells(i_cell) )) = [];
        OFFP(find(OFFP == badcells(i_cell) )) = [];
    end
    
    allcells{i_exp,1}.ONP = ONP ; 
    allcells{i_exp,1}.OFFP = OFFP;
end
%
eval(sprintf('save %s/allcells.mat allcells', savedir))

%% GLM scores
clear;
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat',BD.Cell_Selection))
GLMscores = allcells;


GLMType.debug = false;
GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
GLMType.nullpoint = 'mean'; 
GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
GLMType.CONVEX = true; % with relation to the filters .. are parameters used linearly in the GLM. 
GLMType.specialchange = false;
GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
GLMType.fitname  = GLM_fitname(GLMType)