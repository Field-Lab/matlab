clear datarun

%% load data
opt=struct('verbose',1,'load_params',1,'load_neurons',1);
datarun{1} = load_data('/Volumes/Mahler/Analysis/gfield/2010-03-05-2/data013/data013', opt);
datarun{2} = load_data('/Volumes/Mahler/Analysis/gfield/2010-03-05-2/data016-from-data013/data016-from-data013', opt);

% map cells
datarun=map_cell_types(datarun);

% parse stimulus
datarun{2}.names.stimulus_path = '/Volumes/Mahler/Analysis/gfield/2010-03-05-2/stimuli/s16';
datarun{2}=load_stim(datarun{2}); 


%% Look at raw spike rasters
epoch_condition = 28;
epoch_indices = find(datarun{2}.stimulus.trial_list == epoch_condition);


raster(datarun{2}, datarun{2}.cell_types{1}.cell_ids(4), datarun{2}.stimulus.triggers(epoch_indices), 'stop', 10);

%%
if 1 %CRG f1/f2  screen
   clf; set(gcf, 'color', 'white');

   index=get_cell_indices(datarun{2},{4});

   nr=2;
   k=1; kmin=1; kmax=length(index); ha=loop_slider(k,kmin,kmax);
   while k
      k=round(get(ha,'Value')); 

       r=zeros(2,length(datarun{nr}.stimulus.combinations), datarun{nr}.stimulus.repetitions);      
       for i=1:length(datarun{nr}.stimulus.combinations)

           t1=find(datarun{nr}.stimulus.trial_list==i);
           for ii=1:length(t1)
               bin=.01;
               t=hist(datarun{nr}.spikes{index(k)},[0:bin:8]+datarun{nr}.stimulus.triggers(t1(ii)));
               t=t(2:end-1); 

               nfft=2^nextpow2(length(t)); 
               y=fft(t,nfft)/length(t);
               f=1/bin/2*linspace(0,1,nfft/2);

               [junk t3]=min(abs(f-120/datarun{nr}.stimulus.trials(i).TEMPORAL_PERIOD));
               [junk t4]=min(abs(f-120/datarun{nr}.stimulus.trials(i).TEMPORAL_PERIOD*2));
               r(1,i,ii)=2*abs(y(t3));
               r(2,i,ii)=2*abs(y(t4));  
           end         
       end

       rr=zeros(length(datarun{nr}.stimulus.params.SPATIAL_PERIOD),2,3);      
       for i=1:length(datarun{nr}.stimulus.params.SPATIAL_PERIOD)

           t=find([datarun{nr}.stimulus.combinations.SPATIAL_PERIOD]==datarun{nr}.stimulus.params.SPATIAL_PERIOD(i));

           rr(i,1,:)=max(r(1,t,:));
           rr(i,2,:)=max(r(2,t,:));

       end

       plot((1000/5.8)./fliplr(datarun{nr}.stimulus.params.SPATIAL_PERIOD),flipud(mean(rr(:,1,:),3)),'k')
       hold on
       plot((1000/5.8)./fliplr(datarun{nr}.stimulus.params.SPATIAL_PERIOD),flipud(mean(rr(:,2,:),3)),'r')
       hold off
       set(gca, 'xscale', 'log') 

       xlabel('spatial (cpd)')
       ylabel('amplitude(spike/s)')
       title(sprintf('cell-id=%d',datarun{nr}.cell_ids(index(k)))); 

       uiwait;
   end
end