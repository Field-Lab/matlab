function delete_spikes(thr)
global spikes spikeTimes hist_fig

display('You have removed all spikes above the threshold. Load them again if you want them back')

% which spikes to delete
spikes2delete=min(spikes)>=thr;

% delete spikes
spikes(:,spikes2delete)=[];

% total amount of spikes by stim files
cumSpikes=zeros(1,length(spikeTimes));
cumSpikes(1)=length(spikeTimes{1});
for i=2:length(spikeTimes)
    cumSpikes(i)=cumSpikes(i-1)+length(spikeTimes{i});
end

% delete spike timestamps of invalidated spikes
spikeTimes{1}(spikes2delete(1:cumSpikes(1)))=[];
for i=2:length(spikeTimes)
    tmp=spikes2delete(cumSpikes(i-1)+1:cumSpikes(i));
    spikeTimes{i}(tmp)=[];
end
figure(hist_fig)
deleteSpikesThrUI=uicontrol('style','text','string','removed!','BackgroundColor','m');
set(deleteSpikesThrUI,'Units', 'normalized','Position',[0.4 0.49 0.2 0.04],'fontsize',10,'fontweight','bold');

