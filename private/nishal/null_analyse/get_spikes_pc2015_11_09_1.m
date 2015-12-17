
function [spksGen,spksGen_hr]=get_spikes_pc2015_11_09_1(neuronPath,cellID,stim_length)

neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
CellSpkTimes=double(neuronFile.getSpikeTimes(cellID))/20000;
TTL=double(neuronFile.getTTLTimes())/20000;

spks=double(CellSpkTimes);


icell=0;

%% Load spikes 
      
        % Spike loading
        spikes=spks;
        
        % Align the spikes and the movies;
        spikes_adj=spikes;
        n_block=0;
        for i=1:(length(TTL)-1)
            actual_t_start=TTL(i);
            supposed_t_start=n_block*100/120;
            idx1=spikes > actual_t_start;
            idx2=spikes < TTL(i+1);
            spikes_adj(find(idx2.*idx1))=spikes(find(idx2.*idx1))+supposed_t_start-actual_t_start;
            n_block=n_block+1;
        end
        clear spikes
        spike.home=spikes_adj;
        clear spikes_adj;
        %%
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
        
end