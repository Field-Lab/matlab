
startup_null_analyse_tenessee

% %%
% datarun.names.rrs_neurons_path='/Volumes/Analysis/2012-08-09-3/data002/data002.neurons';
%     
% opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
% datarun=load_data(datarun,opt);
% 
% triggers=datarun.triggers; %onsets of the stimulus presentation
% 
% mdf_file='/Volumes/Analysis/movie-xml/RGB-8-1-0.48-11111.xml';
% 
%  [mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
%     triggers, 1,2);
% 
% [mvi] = load_movie(mdf_file, triggers);
% 
% % figure
% % imagesc(mov(:,:,2))
% 
% cellID=find(datarun.cell_ids==2926)
% 
% 
% spikes=datarun.spikes{cellID};
% 
% % spikes are in s. convert to ms
% spikes=round(spikes*1000);
% %fr - frames - are in ms
% %fr=round(triggers(1)*1000:refresh:triggers(end)*1000);
% % make fr better
% 
% fr=[];
% for itrig=1:length(triggers)
% fr=[fr,triggers(itrig)*1000+refresh*[0:99]];
% end
% fr=fr';
% 
% sta=zeros(height,width,30); %height, width, frames back
% tic
% icnt=0;
% for i=spikes'
%  
%     start=find(fr>i,1)-30; 
%     if(start>1000)
%     icnt=icnt+1
%         for j=1:30
%         F = round(mvi.getFrame(start+j).getBuffer);
%         sta(:,:,j) = sta(:,:,j) + reshape(F(1:3:end),width,height)'+reshape(F(2:3:end),width,height)'+reshape(F(3:3:end),width,height)';
%         end
%     end
% end
% sta=sta/icnt;
% 
% toc  
% 
% % See STA .. STA should be good and its a proof that the code is working
% % fine ..
% 
% figure;
% for j=1:30
% imagesc(sta(:,:,j));
% colormap gray
% caxis([min(sta(:)),max(sta(:))]);
% colorbar
% pause(1/120);
% end


%% 

datafile = '2012-08-09-3/data002';
type_name= cell(1,1);
type_name{1}='On Parasol';

opt=struct('load_all',true);
datarun=load_data(datafile,opt)
datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);
OnParasol_vis_id=datarun.cell_types{1}.cell_ids;

neuronPairsRefVsNew = crossIdentifyNeuronIDs('/Volumes/Analysis/2012-08-09-3/data002', '/Volumes/Analysis/2012-08-09-3/data000',OnParasol_vis_id);
ref_on_parasol_cells=neuronPairsRefVsNew(:,2);

%% Select one cell
vision_id=4052;
ref_cell_idx=find(neuronPairsRefVsNew(:,1)==vision_id);
ref_vision_id=neuronPairsRefVsNew(ref_cell_idx,2);


matlab_ids=find(datarun.cell_ids==vision_id);
sta_data002=datarun.stas.stas{matlab_ids(1)};


%%
clear datarun
% datarun.names.rrs_neurons_path='/Volumes/Analysis/2012-08-09-3/data000/data000.neurons';
%     
% opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true,'load_all',1);
% datarun=load_data(datarun,opt);
datafile = '2012-08-09-3/data000';


opt=struct('load_all',true);
datarun=load_data(datafile,opt)
datarun=load_sta(datarun)
datarun=load_params(datarun)
datarun=load_ei(datarun,'all','array_type',519);

triggers=datarun.triggers; %onsets of the stimulus presentation



%%
% matlab_ids=zeros(length(ref_on_parasol_cells),1);
% for icell=1:length(ref_on_parasol_cells)
% matlab_ids(icell)=find(datarun.cell_ids==ref_on_parasol_cells(icell));
% end

matlab_ids=find(datarun.cell_ids==ref_vision_id);
sta=datarun.stas.stas{matlab_ids(1)};
[sig_stixels, params, rf_strength_out] = significant_stixels(sta)
[sig_row,sig_col]=find(sig_stixels==1);
sig_pix=[sig_row,sig_col];
%%
mdf_file='/Volumes/Analysis/movie-xml/BW-20-5-0.48-11111-30x30-60.35.xml';

 [mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);

% 'refresh' parameter is the amount of time a particular image
% is shown before it is changed.
subdivideRefresh= 1; % So, each frame refresh is subdivided into this many parts. 


% For example , if you are showing BW-10-4 .. The refresh for the movie
% would be 8.33*4 ms .. But we still want to have temporal resolution of
% 8.33ms. So, set subdivideRefresh=4. If you want temporal resolution of
% 4.16ms, then set subdivideRefresh=8. subdivideRefresh must be an integer

distinctImagesbetweenTriggers=round(mean(abs(diff(triggers)))*1000/refresh);
% 'distinctImagesbetweenTriggers' is the number of distinct images which
% are shown between refreshes. 

fr=[];
for itrig=1:length(triggers)
fr=[fr,triggers(itrig)*1000+(refresh/subdivideRefresh)*[0:distinctImagesbetweenTriggers*subdivideRefresh-1]];
end
fr=fr';
% fr are times in ms, when a new bin should start. At times when the movie
% frame changes, a new bin should also start. So, for each bin, there is
% only one movie frame. 

frMov=[];
for itrig=1:length(triggers)
frMov=[frMov,triggers(itrig)*1000+(refresh)*[0:distinctImagesbetweenTriggers-1]];
end
frMov=frMov';


% frMov are times when the frames change. frMov should times
% should be every subdivideRefresh-th of fr times.



interestingCellIds=matlab_ids;%[1:length(datarun.cell_ids)];% These would be matlabcellIds - NOT.. datarun.cell_ids;% This would change in future, I think
nCells=length(interestingCellIds);

spk_coll=cell(nCells,1);
for icell=1:nCells
    icell
    % Next, count spikes in each of these bins.
cellID=interestingCellIds(icell);
spikes=datarun.spikes{cellID};
% spikes are in s. convert to ms
spikes=round(spikes*1000);
spks=zeros(length(fr)-1,1);
for ibin=1:length(fr)-1
spks(ibin)=sum(double(spikes>=fr(ibin) & spikes<fr(ibin+1)));
end
spk_coll{icell}=spks;
end


% Make movie  
mov=zeros(height,width,3,length(fr)-1);
for ibin=1:length(fr)-1
frameNo=find(fr(ibin)>=frMov,1,'last')-1; 
if(frameNo<duration)
 F = round(mvi.getFrame(frameNo).getBuffer);
 mov(:,:,1,ibin)= reshape(F(1:3:end),width,height)';
 mov(:,:,2,ibin)= reshape(F(2:3:end),width,height)';
 mov(:,:,3,ibin)= reshape(F(3:3:end),width,height)';
else
    break
end
end

% Correct for short movie
% If there was an error because movie was small in length, correct for it.
finalLen=ibin;
mov=mov(:,:,:,1:finalLen);
for icell=1:nCells
   spk_coll{icell}=spk_coll{icell}(1:finalLen); 
end

% Play movie ? 
figure
for ibin=1:20
imagesc(mov(:,:,:,ibin));
title(sprintf('%d',ibin));
pause
end

