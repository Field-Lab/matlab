
startup
%%
% I have code for single electrode only .. how to exploit multiple electrodes? 
cell_no=168;
patt = 23;
%%
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile('/Volumes/Analysis/2012-09-24-3/data007/data007.ei');
ids = eiFile.getIDList();

ei = eiFile.getImage(cell_no,0);  % (id, errorType (Standard Deviation of the Mean: 0, Variance of the Mean: 1))
% for standard version of matlab, no second argument (standard deviation vs. variance not an option)

% returns ei(average:1 error:2, electrode + 1, time index)
% 3 dimensional array with first dimension referring to whether value is the average voltage (1) or the error (2), second dimension referring to electrode number, and % third dimension referring to point in time

maxElectrode = eiFile.getMaxElectrode(ei);

%% Which pattern to use ? 
% Get which patterns uses maxElectrode
pathToAnalysisData = '/Volumes/Analysis/2012-09-24-3/data008/';
patt_list=[];
for patternNo=1:368 %patt:patt
% Find movie indices
movieNos = [];
patternNoString = ['p' num2str(patternNo)];
files = dir([pathToAnalysisData patternNoString]);

for i = 1:length(files)
   if strfind(files(i).name, patternNoString) == 1
       mIndices = strfind(files(i).name, 'm');
       movieNos = [movieNos str2double(files(i).name(mIndices(end)+1:end))]; %#ok<AGROW>
   end
end
movieNos = sort(movieNos);


load(sprintf('/Volumes/Analysis/2012-09-24-3/data008/pattern_files/pattern%d_m%d',patternNo,movieNos(1)));
elects = Pattern.channel;
if(ismember(maxElectrode,elects,'rows'))
patt_list=[patt_list;patternNo];
end

end

%% 
maxElectrode = elects(1);
%% Find cells near current electrode 

%% Movie number
patternNo =24;%patt_list(1);

movieNos = [];
patternNoString = ['p' num2str(patternNo)];
files = dir([pathToAnalysisData patternNoString]);

for i = 1:length(files)
   if strfind(files(i).name, patternNoString) == 1
       mIndices = strfind(files(i).name, 'm');
       movieNos = [movieNos str2double(files(i).name(mIndices(end)+1:end))]; %#ok<AGROW>
   end
end
movieNos = sort(movieNos);

movieNo = movieNos(2);

%% 
no_cells=1;
wave_sample_len=size(ei,3);
waveforms=zeros(wave_sample_len,no_cells+1);
waveforms(:,2)=reshape(ei(1,maxElectrode+1,:),[1,size(ei,3)]);

% Make snippets
snippets=struct;
no_snippets=1;
snippet_sz=10;
snippets(no_snippets).data=[];

% 
icount=1;

for isz=1:snippet_sz:no_snippets*snippet_sz
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(sprintf('/Volumes/Analysis/2012-09-24-3/data008/p%d',patternNo),[],0,patternNo,movieNo,20000);

data=[];
for iisz=1:snippet_sz
data =[data ;[reshape(DataTraces(iisz,maxElectrode,:),[100,1]);reshape(DataTraces(iisz,maxElectrode,1),[1 1])*ones(100,1)]];
end


idx=[1:200*snippet_sz];
data_event=0*data;
data_event(mod(idx,200)==1) =1; 
data=data';
data_event=data_event';

data=data(:,1:end-1);
data=data-mean(data);
data_event=data_event(:,1:end-1);



snippets(icount).data=data;
snippets(icount).data_event=data_event;

icount=icount+1;
end
time_sample_len=200*snippet_sz-1;

figure;
subplot(3,1,1)
plot(waveforms);
subplot(3,1,2)
plot(reshape(DataTraces(:,maxElectrode,:),[15,100])')
subplot(3,1,3)
plot(snippets(1).data)


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

%% CVX _polar - ADMM
%debug
no_snippets=no_snippets;

lam=5.5%5.5;%0.5

stim_est=zeros(wave_sample_len,no_snippets);
y = stim_est;
mu=0.1; %0.001 ?? - lower is better ?? 0.01
stim_avg=zeros(wave_sample_len,1);

% initialize stim_avg using clustering
% offset business might have some problems??
eps=0.001; % worked for 0.001, 0.0001
gam=100;
for dist_iter=1:1
    cvx_optval=[];
    for  snippet_select=1:no_snippets
        
        if(dist_iter==1)
            snippets(snippet_select).opt_log=[];
            snippets(snippet_select).stim_est_log=[];
        end
        
        reweight_exp=1.5*ones(vec_len,1);
        lam=reweight_exp;
         
        A_stim = sparse(time_sample_len,wave_sample_len);
            for i_time=1:time_sample_len
                if(snippets(snippet_select).data_event(i_time)==1)
                    A_stim(i_time:i_time+wave_sample_len-1,1:end) = A_stim(i_time:i_time+wave_sample_len-1,1:end) + eye(wave_sample_len);
                end 
            end
        snippets(snippet_select).A_stim=A_stim;
        
        for i_reweight=1:5
            [dist_iter snippet_select i_reweight]
            % Another loop for iterated L1?
            
            
           
            
            [C,U,V,stim_est(:,snippet_select),cvx_optval(snippet_select)]=cvx_polar_admm(snippets(snippet_select).data,snippets(snippet_select).data_event,radii_vec,theta_vec,time_sample_len,dict,lam,vec_len,wave_sample_len,stim_avg,mu,y(:,snippet_select),gam,A_stim);
            
            lam=reweight_exp./(eps+abs(C));
            snippet_select
            % end`
           
        end
       
    snippets(snippet_select).C=C;
    snippets(snippet_select).U=U;
    snippets(snippet_select).V=V;
    snippets(snippet_select).opt_log=[snippets(snippet_select).opt_log;cvx_optval(snippet_select)];
    snippets(snippet_select).stim_est_log=[snippets(snippet_select).stim_est_log,stim_est(:,snippet_select)];

    
    end
    
    stim_avg=mean(stim_est,2) + mu*mean(y,2);
    y=y+(1/mu)*(stim_est - repmat(stim_avg,1,no_snippets));
    
  
    
    figure;
    plot(stim_avg)
    title(sprintf('Iteration %d',dist_iter));
    hold on
    plot(waveforms(:,1),'r')
    title('Artifact and it estimation')
    legend('Estimated artifact','Artifact')
    hold on
    plot(stim_est)
          % pause
            
end

%%

subplot(2*no_snippets+1,1,1);
stim_avg=mean(stim_est,2);
plot(stim_avg)



    for isnip=1:no_snippets
        
        
snip_select=isnip;

data=snippets(snip_select).data;
data_event=snippets(snip_select).data_event;
%cells_fire=snippets(snip_select).cells_fire;
A_stim = snippets(snip_select).A_stim;


        subplot(2*no_snippets+1,1,2*(isnip-1)+2);
plot(data)
title('Data - 0 mean');
       
        
    subplot(2*no_snippets+1,1,2*(isnip-1)+3);
    
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


    
    
  %  plot(recons);
    
    y=data_event*0;
    x=data_event*0;
    x(round(spike_times{icell} + 0.5))=spike_amps{icell};
    y(round(spike_times{icell}+0.5))=1;
    hold on;
    stem(y.*x,'g');
    title(sprintf('Estimated spike train , Snippet: %d',isnip));
    legend('recons','CBP')
    
    end
    
    

