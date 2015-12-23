

datarun.names.rrs_neurons_path='/Volumes/Analysis/2012-08-09-3/data002/data002.neurons';
    
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

mdf_file='/Volumes/Analysis/movie-xml/RGB-8-1-0.48-11111.xml';

 [mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);

% figure
% imagesc(mov(:,:,2))

cellID=find(datarun.cell_ids==2926)


spikes=datarun.spikes{cellID};

% spikes are in s. convert to ms
spikes=round(spikes*1000);
%fr - frames - are in ms
%fr=round(triggers(1)*1000:refresh:triggers(end)*1000);
% make fr better

fr=[];
for itrig=1:length(triggers)
fr=[fr,triggers(itrig)*1000+refresh*[0:99]];
end
fr=fr';

sta=zeros(height,width,30); %height, width, frames back
tic
icnt=0;
for i=spikes'
 
    start=find(fr>i,1)-30; 
    if(start>1000)
    icnt=icnt+1
        for j=1:30
        F = round(mvi.getFrame(start+j).getBuffer);
        sta(:,:,j) = sta(:,:,j) + reshape(F(1:3:end),width,height)'+reshape(F(2:3:end),width,height)'+reshape(F(3:3:end),width,height)';
        end
    end
end
sta=sta/icnt;

toc  