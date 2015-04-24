clear
datapath='/Volumes/Lab/Users/Nora/NSEM/GLM_Output/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/';'2012-09-27-3/';'2013-08-19-6/';'2013-10-10-0/'];
fit_type{1}='WN';
fit_type{2}='NSEM';
fittypepath{2}=[fit_type{2} '_mapPRJ/'];
fittypepath{1}=[fit_type{1} '_mapPRJ/'];

for fittype=1
    for exp=4
        
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

%%
K_time_avg = mean(K_time,1);
K_time_avg_ON = mean(K_time(type == 1 , :),1);
K_time_avg_OFF = mean(K_time(type == 2 , :),1);
PS_avg_ON = mean(PS(type == 1 , :),1);
PS_avg_OFF = mean(PS(type == 2 , :),1);
mu_avg_ON = mean(mu(type == 1));
mu_avg_OFF = mean(mu(type == 2));

%%
% experiment dependent
[SPars, ~, ~, ~] = Directories_Params_v23(exp_names(exp,1:end-1), fit_type{fittype}, 'mapPRJ');
[testmovie] = loadmoviematfile(exp_names(exp,1:end-1), fit_type{fittype}, '8pix_Identity_8pix','testmovie');

%%
% replace common filters
for file=1:n_cells
    disp(file)
    load([datapath fittypepath{fittype} exp_names(exp,:) matfiles(file).name]);
    if type(file) == 1
        fittedGLM.linearfilters.Stimulus.time_rk1 = K_time_avg_ON';
        % fittedGLM.linearfilters.TonicDrive.Filter = mu_avg_ON';
        % fittedGLM.linearfilters.PostSpike.Filter = PS_avg_ON';
    else
        fittedGLM.linearfilters.Stimulus.time_rk1 = K_time_avg_OFF';
        % fittedGLM.linearfilters.TonicDrive.Filter = mu_avg_OFF';
        % fittedGLM.linearfilters.PostSpike.Filter = PS_avg_OFF';
    end
    fittedGLM.linearfilters.Stimulus.Filter = reshape(fittedGLM.linearfilters.Stimulus.space_rk1(:)*fittedGLM.linearfilters.Stimulus.time_rk1',13,13,30);
    fittedGLM.GLMType.Subunits = false;
    load(['/Volumes/Lab/Users/Nora/NSEM/BlockedSpikes/' exp_names(exp,:) fit_type{fittype} '_mapPRJ/organizedspikes_' fittedGLM.cellinfo.cell_savename '.mat'], 'organizedspikes')
    xval=eval_xvalperformance_NEW_CP(fittedGLM, SPars.slv, organizedspikes, 0, testmovie);
    bps_new(file)=xval.glm_normedbits;
end

%%
plot(bps, bps_new, 'o')
hold on
plot([0 1], [0 1])
xlim([0,0.75])
ylim([0,0.75])
axis square
xlabel('BPS with Individual Filters')
ylabel('BPS with Average Filters')
title('Average K')



