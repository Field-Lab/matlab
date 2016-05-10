function [linear_drive, res] = linear_filter_drive(movie)
Conv = 20;
Cutoff = -0.1; % BPS cutoff
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_noCP_p8IDp8/standardparams/NSEM_mapPRJ/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
n_cells = 40;
res.cells = zeros(n_cells, 1);
res.centers = zeros(n_cells, 2);
i_cell = 0;

for exp = 1:2
    matfiles = dir([datapath exp_names(exp,:) 'conv_blocks_' num2str(Conv) '/*.mat']);
        
    for file = 1:length(matfiles)
        
        % Load FittedGLM
        load([datapath exp_names(exp,:) 'conv_blocks_' num2str(Conv) '/' matfiles(file).name]);
        
        if fittedGLM.xvalperformance.glm_normedbits > Cutoff  % Only use the fit if it has relatively converged
            i_cell = i_cell + 1
            
            % Save cell info
            res.cells(i_cell) = fittedGLM.cellinfo.cid;
            res.centers(i_cell,1) = fittedGLM.cellinfo.slave_centercoord.y_coord;
            res.centers(i_cell,2) = fittedGLM.cellinfo.slave_centercoord.x_coord;
            
            % Load a little info
            fittedGLM.center_coord = fittedGLM.cellinfo.slave_centercoord;
            fittedGLM.inputstats.mu_avgIperpix = 64;
            fittedGLM.inputstats.range = 255;
          
            
            % Get the linear drive: stimulus + tonic
            for block = 1:length(movie)
                linear_drive{i_cell,block} = glm_linear_predict(fittedGLM, movie{block}.matrix);
            end
        end
        
    end
    
end
end