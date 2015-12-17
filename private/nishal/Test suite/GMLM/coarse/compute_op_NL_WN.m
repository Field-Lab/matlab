function fitGMLM_log = compute_op_NL_WN(WN_datafile,movie_xml,stim_length,cellIDs,user_STA_depth,destination,contrast_factor,save_location)


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);
dt=1/120;
%% Find stimulus and response
% cellID=1531;
for cellID = cellIDs
% user_STA_depth=30;
cellID
%extract_movie_response2;

   master_idx         = find(datarun.cell_ids == cellID);
    
      
        % Spike loading
        spikes=datarun.spikes{master_idx};
    
        % make STA 3D ? 
        %glm_cellinfo.WN_STA = squeeze(sum(glm_cellinfo.WN_STA,3)); % Doubt!!!!!!!
        clear cell_savename
        
        % Align the spikes and the movies;
        spikes_adj=spikes;
        n_block=0;
        for i=1:(length(datarun.triggers)-1)
            actual_t_start=datarun.triggers(i);
            supposed_t_start=n_block*100/120;
            idx1=spikes > actual_t_start;
            idx2=spikes < datarun.triggers(i+1);
            spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
            n_block=n_block+1;
        end
        clear spikes
        spike.home=spikes_adj;
        clear spikes_adj;
        %
        spksGen = zeros(stim_length*120,1);
        for ispike=1:length(spike.home)
            spksGen(floor(spike.home(ispike)*120)+1)=1;
        end
        spksGen = spksGen(1:stim_length*120);
        
        spksGen_hr = zeros(stim_length*1200,1);
        for ispike=1:length(spike.home)
            spksGen_hr(floor(spike.home(ispike)*1200)+1)=1;
        end
        spksGen_hr = spksGen_hr(1:stim_length*1200);
%% Load fits

ASM_link =['/Volumes/Lab/Users/bhaishahster/',destination,sprintf('/Cell_%d.mat',cellID)];
data1 = load(ASM_link);
fitGMLM_log = data1.fitGMLM_log;


%%
close all
clear fit_nl
for nSU=1:length(fitGMLM_log)
    nSU
fitGMLM = fitGMLM_log{nSU};
% 
% % compute NL
% [g,meanl,meanR] = compute_op_NL(fitGMLM,maskedMovdd,spksGen_hr,dt);
% fitGMLM_log{nSU}.NL_op = g;
% fitGMLM_log{nSU}.NL.meanl = meanl;
% fitGMLM_log{nSU}.NL.meanR = meanR;
% fitGMLM_log{nSU}.NL.dt = dt;
nbins=10;
lam = fitGMLM.data_act.lam;
[fit_params, gof,fitres1,fitres] = fit_nonlinearity_alex( spksGen/dt, lam,nbins)
fit_nl(nSU).fit=fitres;
fit_nl(nSU).fit1 = fitres1;
fit_nl(nSU).fit1_params = fit_params;
end

data1.fitGMLM_log = fitGMLM_log;

save_folder =['/Volumes/Lab/Users/bhaishahster/',save_location,sprintf('/Cell_%d.mat',cellID)];
save(save_folder,'fit_nl','-v7.3');
end

end