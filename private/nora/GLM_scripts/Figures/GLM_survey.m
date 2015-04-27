clear
datapath='/Volumes/Analysis/nora/NSEM/GLM_Output/rk1_MU_PS_CP_p8IDp8_orig/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

for fittype=2
    for exp=1
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_cells=length(matfiles);
        
        % Collect info from files
        for file=1:n_cells
            disp(file)
            load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
            % rasters{file} = fittedGLM.xvalperformance.rasters;
            plotraster(fittedGLM.xvalperformance, fittedGLM);
            pause()
        end
    end
end


