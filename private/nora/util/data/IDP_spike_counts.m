function spike_counts = IDP_spike_counts(PSTH, image_transitions)
% PSTH and image transitions should be in the same units, i.e. frames or
% time or bins or whatever
n_images = length(image_transitions);
spike_counts = zeros(n_images,1);
for i = 1:n_images
    try
        temp = PSTH(image_transitions(i):image_transitions(i+1));
    catch
        temp = PSTH(image_transitions(i):end);
    end
   spike_counts(i) = sum(temp);
end

end