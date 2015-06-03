clear
datapath='/Volumes/Analysis/nora/NSEM/GLM_Output/old_fits/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fittypepath{2}='NSEM_mapPRJ/';
fittypepath{1}='WN_mapPRJ/';

for fittype=1
    for exp=1
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(exp,:) '*.mat']);
        n_cells=length(matfiles);
        
        %Initialize matrices
        mu=zeros(n_cells,1);
        bps=zeros(n_cells,1);
        type=zeros(n_cells,1);
        PS=zeros(n_cells,120);
        K_space=zeros(n_cells, 13, 13);
        K_time=zeros(n_cells, 30);
        
        % Collect info from files
        for file=1:n_cells
            disp(file)
            load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
            K_space(file,:,:)=fittedGLM.linearfilters.Stimulus.space_rk1;
            K_time(file,:)=fittedGLM.linearfilters.Stimulus.time_rk1;
            mu(file)=fittedGLM.linearfilters.TonicDrive.Filter;
            PS(file,:)=fittedGLM.linearfilters.PostSpike.Filter;
            bps(file)=fittedGLM.xvalperformance.glm_normedbits;
            if strcmp(matfiles(file).name(2), 'N')
                type(file)=1;
            else
                type(file)=2;
            end
        end
    end
end
