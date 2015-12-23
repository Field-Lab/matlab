function []=fcn_real_anal(istim,movie_idx,movieNo, patternNo,cell_no, ei,maxElectrode,mu,lam_re,eps,gam);


%% 
snip_len=138;

no_cells=1;
data_wave_len=size(ei,3);
wave_sample_len=snip_len-1;%min(snip_len-1,data_wave_len);
waveforms=zeros(wave_sample_len,no_cells+1);
waveforms(1:data_wave_len,2)=reshape(ei(1,maxElectrode+1,1:data_wave_len),[1,data_wave_len]);
% natural error if cell spike waveform length > wave_sample_len
% Make snippets
snippets=struct;
no_snippets=1;
snippet_sz=15;
snippets(no_snippets).data=[];

% 
icount=1;

for isz=1:snippet_sz:no_snippets*snippet_sz
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(sprintf('/Volumes/Analysis/nishal/data008/p%d',patternNo),[],0,patternNo,movieNo,20000);

data=[];
for iisz=1:snippet_sz
%data =[data ;[reshape(DataTraces(iisz,maxElectrode,:),[100,1]);reshape(DataTraces(iisz,maxElectrode,1),[1 1])*ones(100,1)]];
data =[data ;[reshape(DataTraces(iisz,maxElectrode,1:snip_len),[snip_len,1])]];

end


idx=[1:snip_len*snippet_sz];
data_event=0*data;
data_event(mod(idx,snip_len)==1) =1; 
data=data';
data_event=data_event';

data=data(:,1:end-1);
mdata=mean(data);
data=data-mean(data);
data_event=data_event(:,1:end-1);



snippets(icount).data=data;
snippets(icount).data_event=data_event;

icount=icount+1;
end
time_sample_len=snip_len*snippet_sz-1;

fig2=figure;
subplot(3,1,1)
plot(waveforms);
subplot(3,1,2)
plot(reshape(DataTraces(1:15,maxElectrode,:),[15,200])')
subplot(3,1,3)
plot(snippets(1).data)

% Load elecResp
load(sprintf('/Volumes/Analysis/2012-09-24-3/data008/elecResp_n%d_p%d.mat',cell_no,patternNo));
lats=elecResp.analysis.latencies{movie_idx};
true_resp=[];
for ilat=1:snippet_sz % for 15 trials
    datx=zeros(1,snip_len);
    if(lats(ilat)~=0)
    datx(round(lats(ilat)))=1;
    end
    true_resp=[true_resp,datx];
end

%% Artifact estimate
dat_art_xx= elecResp.analysis.estArtifact{movie_idx};
waveforms(1:length(dat_art_xx),1) = dat_art_xx-mdata;
%%
% waveforms as cells to be used to make dictionary
init_features = cell(no_cells,1);
for icell=1:no_cells
    init_features{icell} = waveforms(:,icell+1);
end


% Code parts to do CBP
accuracy=0.1;
spacings = ...
    polar_1D_delta(init_features, accuracy);

[radii, thetas] = ...
    polar_1D_get_radii_angles(init_features, spacings);
snippet_select=1;
data_cell=cell(1);
data_cell{1}=snippets(snippet_select).data';

[gps dicts] = ...
    PrecomputeDictionaries(init_features, spacings, data_cell);


dict=dicts{time_sample_len};

vec_len =size(dict,2)/3;



num_blocks_per_feature = cellfun(@(C) size(C,1), gps{time_sample_len});
%lambda_vec = multirep(lambda(:), num_blocks_per_feature);
radii_vec = multirep(radii(:), num_blocks_per_feature);
theta_vec = multirep(thetas(:), num_blocks_per_feature);


%%
%% CVX _polar - ADMM
%debug
no_snippets=no_snippets;

%lam_re=0.25;%0.5

stim_est=zeros(wave_sample_len,no_snippets);
y = stim_est;
%mu=100; % ?? - lower is better ?? 0.01
stim_avg=zeros(wave_sample_len,1);
%stim_avg=artifact_mean(1:wave_sample_len);
% initialize stim_avg using clustering
% offset business might have some problems??
%eps=0.1; % worked for 0.001, 0.0001
%gam=1%1;
for dist_iter=1:1
    cvx_optval=[];
    for  snippet_select=1:no_snippets
        
        if(dist_iter==1)
            snippets(snippet_select).opt_log=[];
            snippets(snippet_select).stim_est_log=[];
        end
        
        reweight_exp=lam_re*ones(vec_len,1);
        lam=reweight_exp;
         
        A_stim = sparse(time_sample_len,wave_sample_len);
            for i_time=1:time_sample_len
                if(snippets(snippet_select).data_event(i_time)==1)
                    A_stim(i_time:i_time+wave_sample_len-1,1:end) = A_stim(i_time:i_time+wave_sample_len-1,1:end) + eye(wave_sample_len);
                end 
            end
        snippets(snippet_select).A_stim=A_stim;
        
        for i_reweight=1:8
            [movieNo dist_iter snippet_select i_reweight]
            % Another loop for iterated L1?
            
            
           
            
            [C,U,V,stim_est(:,snippet_select),cvx_optval(snippet_select)]=cvx_polar_admm_sqrt2(snippets(snippet_select).data,snippets(snippet_select).data_event,radii_vec,theta_vec,time_sample_len,dict,lam,vec_len,wave_sample_len,stim_avg,mu,y(:,snippet_select),gam,A_stim);
            
            lam=reweight_exp./(eps+abs(C));
            
            % end
            sum(C>0.01)
            sum(C)
           
        end
       
    snippets(snippet_select).C=C;
    snippets(snippet_select).U=U;
    snippets(snippet_select).V=V;
    snippets(snippet_select).opt_log=[snippets(snippet_select).opt_log;cvx_optval(snippet_select)];
    snippets(snippet_select).stim_est_log=[snippets(snippet_select).stim_est_log,stim_est(:,snippet_select)];

    
    end
    
    stim_avg=mean(stim_est,2) + mu*mean(y,2);
    y=y+(1/mu)*(stim_est - repmat(stim_avg,1,no_snippets));
    
  
    
    fig0=figure;
    plot(stim_avg)
    title(sprintf('Iteration %d',dist_iter));
    hold on
    plot(waveforms(:,1),'r')
    title('Artifact and it estimation')
    legend('Estimated artifact','Artifact')
    hold on
    plot(stim_est)
    ylim([-500,500]);
          % pause
            
end


%%  Plot data
fig1 = figure;


subplot(2*no_snippets+1+1,1,1);
stim_avg=mean(stim_est,2);
plot(stim_avg)
hold on
plot(waveforms(:,1),'r')
ylim([-500,500])
title('Artifact and it estimation')
legend('Estimated artifact','Artifact')


    for isnip=1:no_snippets
        
        
snip_select=isnip;

data=snippets(snip_select).data;
data_event=snippets(snip_select).data_event;
%cells_fire=snippets(snip_select).cells_fire;
A_stim = snippets(snip_select).A_stim;


        subplot(2*no_snippets+1+1,1,2*(isnip-1)+2);
plot(data)
title('Data - 0 mean');
       
        
    subplot(2*no_snippets+1+1,1,2*(isnip-1)+3);
    %cells_fire = snippets(isnip).cells_fire; 

    y=data_event*0;
    %y(round(cells_fire{icell}))=1;
    hold on;
    stem(true_resp,'r');
    

    
    
    



%
C=snippets(snip_select).C;
U=snippets(snip_select).U;
V=snippets(snip_select).V;


%raw_coeffs = reshape([C,U,V], [], 1);
raw_coeffs = vec([C, U, V]');
%raw_coeffs=[C;U;V];
% Extract transformation params, magnitudes, raw coefficients.
magnitude_threshold = 0.01;
[transform_params, magnitudes,recons] = ...
    polar_1D_extract_fn_cvx(raw_coeffs, ...
    num_blocks_per_feature, ...
    dict, ...
    radii_vec, ...
    gps{time_sample_len}, ...
    spacings, ...
    thetas, ...
    magnitude_threshold); % last entry is magnitude

spike_amps_cell=cell(1);
spike_times_cell=cell(1);
spike_times_cell{1} = transform_params;
spike_amps_cell{1} = magnitudes;
snip_centers=round(time_sample_len/2);

[spike_times, spike_amps] = ConvertSpikeTimesFromCell(spike_times_cell, ...
    spike_amps_cell, snip_centers);


    
    
    
    
    y=data_event*0;
    x=data_event*0;
    x(round(spike_times{icell} + 0.5))=spike_amps{icell};
    y(round(spike_times{icell}+0.5))=1;
    hold on;
    stem(y.*x,'g');
    xlim([1 size(data,2)]);
    
    hold on;
    indx = [1:size(data,2)];
    stim_indx= indx(data_event==1);
    data_win=0*data_event;
    for iisz=1:snippet_sz
    data_win(stim_indx(iisz):stim_indx(iisz)+100)=1;
    end
    stem(y.*x.*data_win,'k');
   
    title(sprintf('Estimated spike train , Snippet: %d',isnip));
    legend('Human est','CBP rest','CBP wind','location','NorthEastoutside');
    
    end
    
     subplot(2*no_snippets+1+1,1,2*no_snippets+1+1);
     
plot(A_stim*stim_avg,'b')%dict*vec([C, U, V]')) %+A_stim*stim_avg);
hold on
plot(data,'r');
hold on 
plot(dict*vec([C, U, V]'),'g')

title('waveform and its estimation');
legend('artifact reconstruction','measured','cell activity');



%     
     print(fig2,sprintf('/Volumes/Analysis/nishal/2012-09-24-3/stim%d/data_mov%d_lam%f_eps%f_gam%f_mu%f.eps',istim,movieNo,lam_re,eps,gam,mu),'-depsc');
     print(fig0,sprintf('/Volumes/Analysis/nishal/2012-09-24-3/stim%d/art_mov%d_lam%f_eps%f_gam%f_mu%f.eps',istim,movieNo,lam_re,eps,gam,mu),'-depsc');
     print(fig1,sprintf('/Volumes/Analysis/nishal/2012-09-24-3/stim%d/dat_mov%d_lam%f_eps%f_gam%f_mu%f.eps',istim,movieNo,lam_re,eps,gam,mu),'-depsc');
     save(sprintf('/Volumes/Analysis/nishal/2012-09-24-3/stim%d/mov%d.mat',istim,movieNo));
%     

end
