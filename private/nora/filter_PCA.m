
clear

%%
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fit_type{1}='WN';
fit_type{2}='NSEM';
fittypepath{2}=[fit_type{2} '_mapPRJ/'];
fittypepath{1}=[fit_type{1} '_mapPRJ/'];

for Iexp = 4
    for fittype = 1
        count = 0;
        
        % Get file list
        matfiles=dir([datapath fittypepath{fittype} exp_names(Iexp,:) 'ON*.mat']);
        n_cells=length(matfiles);
        
        %Initialize matrices
        mu=zeros(n_cells,1);
        bps=zeros(n_cells,1);
        type=zeros(n_cells,1);
        PS=zeros(n_cells,120);
        K_space=zeros(n_cells, 13, 13);
        K_time=zeros(n_cells, 30);
        CP=zeros(n_cells*6,120);
        
        % Collect info from files
        for file=1:n_cells
            disp(file)
            load([datapath fittypepath{fittype} exp_names(Iexp,:) matfiles(file).name]);
            K_space(file,:,:)=fittedGLM.linearfilters.Stimulus.space_rk1;
            K_time(file,:)=fittedGLM.linearfilters.Stimulus.time_rk1;
            mu(file)=fittedGLM.linearfilters.TonicDrive.Filter;
            PS(file,:)=fittedGLM.linearfilters.PostSpike.Filter;
            bps(file)=fittedGLM.xvalperformance.glm_normedbits;
            for pair = 1:6
                cp_max = max(fittedGLM.linearfilters.Coupling.Filter{pair});
                if  cp_max > 0.5
                    count = count + 1;
                    CP(count, :) = fittedGLM.linearfilters.Coupling.Filter{pair}/cp_max;
                end
            end
        end
    end
end

%%
%% PCA
[coeff, score, ~] = pca(PS);
PS_waveform = coeff(:,1);
init = score(:,1);
plot(-PS_waveform)

