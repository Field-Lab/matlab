function [spikesconvbasis]  = prep_convolvespikes_basis(binned_spikes,basis,bins)
% AKHeitman 2014-05-04
% Parameter independent!
% basis should be a vector of [basis vectors , bins] already binned
% t_bin is used to put spike times into their appropriate bin 
% t_bin is duration of each bin in msecs
% bins  is the maximum bin number


vectors = size(basis,1); vectorbins= size(basis,2);
%offset by 1,  so PS filter starts in the next bin!!! 
binned_spikes = binned_spikes(find(binned_spikes < bins) ) + 1;
max(binned_spikes)

convolvedspikes_base                  = zeros(1,bins);
convolvedspikes_base(binned_spikes) = 1; 

convolvedspikes = zeros(vectors, bins + vectorbins - 1);
for i_vec = 1:vectors
    convolvedspikes(i_vec, :) = conv(convolvedspikes_base, basis(i_vec,:), 'full');
end
convolvedspikes = convolvedspikes(:, 1:bins);    

spikesconvbasis = convolvedspikes;
end





