
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
K_time_avg = mean(K_time,1);
PS_avg = mean(PS,1);
mu_avg = mean(mu);
K_space_avg = mean(K_space, 1);
CP_avg = mean(CP(1:count,:));

%%
% experiment dependent
[SPars, ~, ~, ~] = Directories_Params_v23(exp_names(Iexp,1:end-1), fit_type{fittype}, 'mapPRJ');
[~, inputstats, ~] = loadmoviematfile(exp_names(Iexp,1:(end-1)) , fit_type{fittype}, '8pix_Identity_8pix','fitmovie');
[testmovie0]          = loadmoviematfile(exp_names(Iexp,1:(end-1)) , fit_type{fittype}, '8pix_Identity_8pix','testmovie');
testmovie             = testmovie0{1}.matrix(:,:,SPars.slv.testframes);
%%
% replace common filters
for file=1:n_cells
    disp(file)
    load([datapath fittypepath{fittype} exp_names(Iexp,:) matfiles(file).name]);
    fittedGLM.linearfilters.Stimulus.space_rk1 = K_space_avg;
    fittedGLM.linearfilters.Stimulus.time_rk1 = K_time_avg';
    fittedGLM.linearfilters.TonicDrive.Filter = mu_avg';
    fittedGLM.linearfilters.PostSpike.Filter = PS_avg';
    for i = 1:pair
        fittedGLM.linearfilters.Coupling.Filter{i} = max(fittedGLM.linearfilters.Coupling.Filter{i}) * CP_avg';
    end
    fittedGLM.linearfilters.Stimulus.Filter = reshape(fittedGLM.linearfilters.Stimulus.space_rk1(:)*fittedGLM.linearfilters.Stimulus.time_rk1',13,13,30);
    fittedGLM.GLMType.Subunits = false;
    for i_pair = 1:6
        load(['/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/' exp_names(Iexp,:) fit_type{fittype} '_mapPRJ/organizedspikes_' fittedGLM.cellinfo.pair_savename{i_pair} '.mat'], 'organizedspikes')
        neighborspikes.test{i_pair} = subR_raster(organizedspikes.block, SPars.slv);
    end
    load(['/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/' exp_names(Iexp,:) fit_type{fittype} '_mapPRJ/organizedspikes_' fittedGLM.cellinfo.cell_savename '.mat'], 'organizedspikes')
    testspikes_raster.home = subR_raster(organizedspikes.block, SPars.slv);
    xval=eval_xvalperformance(fittedGLM, testspikes_raster, testmovie, inputstats, neighborspikes.test);
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
title('Average Everything')



