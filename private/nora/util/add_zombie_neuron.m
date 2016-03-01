EL = uint32(197);
cluster = uint32(5);
ID  = uint32(2945);
neurons_file = '/Volumes/Analysis/2016-02-17-8/mVision/data000/data000.neurons';

% check is cluster number and ID are consistent
if cluster ~= mod(ID, 15)
   error('check cluster and ID number') 
end

% load up spike times
load([neurons_file '.mat'])
if log(nSamples/log(2)) > 32
    error('uint32 will not be big enough')
end
spike_times = uint32(neuronSpikeTimes{neuronEls == EL+1 & neuronClusters == cluster});
n_spikes = uint32(length(spike_times));

% rewrite neurons file
nrf = edu.ucsc.neurobiology.vision.io.NeuronFile(neurons_file);
nrf.addNeuron(EL, ID, spike_times, n_spikes);
nrf.close;