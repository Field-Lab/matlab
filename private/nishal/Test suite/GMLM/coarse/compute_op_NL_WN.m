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
for nSU=1:10%length(fitGMLM_log)
    nSU
fitGMLM = fitGMLM_log{nSU};
% 
% % compute NL
% [g,meanl,meanR] = compute_op_NL(fitGMLM,maskedMovdd,spksGen_hr,dt);
% fitGMLM_log{nSU}.NL_op = g;
% fitGMLM_log{nSU}.NL.meanl = meanl;
% fitGMLM_log{nSU}.NL.meanR = meanR;
% fitGMLM_log{nSU}.NL.dt = dt;

nbins=20;
lam = fitGMLM.data_act.lam;

% poisson fit
fit_fcn =fit_non_linearity_poissonll(spksGen/dt, lam,nbins);
fit_nl(nSU).fit_nll = fit_fcn;

% least squares fit
hold on;
[fit_params, gof,fitres1,fitres] = fit_nonlinearity_alex( spksGen/dt, lam,nbins) 
fit_nl(nSU).fit_leastSq=fitres;
fit_nl(nSU).fit1_leastSq = fitres1;
fit_nl(nSU).fit1_params_leastSq = fit_params;
% 
title(sprintf('#SU: %d',nSU));

end

data1.fitGMLM_log = fitGMLM_log;

save_folder =['/Volumes/Lab/Users/bhaishahster/',save_location,sprintf('/Cell_%d.mat',cellID)];
save(save_folder,'fit_nl','-v7.3');
end

end


function fit_fcn =fit_non_linearity_poissonll(asr, gs, nbins)

tmp = sort(gs);
 %gs_bins = [tmp(1:bin_size:end) tmp(end)*1.001];
 
% use this if want to fit for saturation!
gs_bins = [tmp(1):(tmp(end)*1.001-tmp(1))/nbins :tmp(end)*1.001];

% gs_bins = [];
% for ibin=0:nbins
% gs_bins = [gs_bins,prctile(gs,100*ibin/nbins)];
% end

nonlinearity=zeros(length(gs_bins)-1,1);
mean_gs_bin=zeros(length(gs_bins)-1,1);
for k=1:length(gs_bins)-1
    nonlinearity(k)=mean(asr(gs>=gs_bins(k) & gs<gs_bins(k+1)));
    mean_gs_bin(k) = mean(gs(gs>=gs_bins(k) & gs<gs_bins(k+1)));
end
mean_gs_bin = mean_gs_bin(~isnan(nonlinearity));
nonlinearity = nonlinearity(~isnan(nonlinearity));


optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter-detailed',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','off');

 [x,fval,exitflag,output,grad,hessian]  = fminunc(@(x)nll(x,gs',asr),[100,100],optim_struct);

a=x(1);b=x(2);
fit_fcn = @(x) (a*x./(b+x));
  
figure;
plot(mean_gs_bin/120,nonlinearity/120,'b--*');
hold on;
plot(mean_gs_bin/120,fit_fcn(mean_gs_bin)/120,'g--','LineWidth',2);

end

function op = rate(inp,params)
a = params(1);
b = params(2);

op = a*inp./(b+inp);
end

function [nllop,grad] = nll(params,inp,op)
op=op/120;
r = rate(inp,params) ;
nllop = (sum(r/120) - op'*log(r))/length(op);

grad = gradnll(params,inp,op);

end

function grad = gradnll(params,inp,op)
r = rate(inp,params);
a = params(1);b=params(2);

grad=0*params;

grad(1) = sum((1/a)*(r.*((1/120) - op./r)));
grad(2) = sum(-r./(b+inp) .*((1/120)  - op./r));

end
