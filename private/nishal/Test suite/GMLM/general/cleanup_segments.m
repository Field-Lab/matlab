function [mov_fit2,binnedSpikeResponses_tr_use]=cleanup_segments(mov_fit,binnedSpikeResponses,cut_time,nSeg)

fr_per_seg = size(mov_fit,2)/nSeg;
new_fr_per_seg = (fr_per_seg - cut_time);
mov_fit2 = zeros(size(mov_fit,1),new_fr_per_seg);

for imov= 1:nSeg
    mov_fit2(:,(imov-1)*new_fr_per_seg+1:imov*new_fr_per_seg) = mov_fit(:,(imov-1)*fr_per_seg + cut_time+1 :imov*fr_per_seg);
    binnedSpikeResponses_tr_use((imov-1)*new_fr_per_seg+1:imov*new_fr_per_seg) = binnedSpikeResponses((imov-1)*fr_per_seg + cut_time+1 :imov*fr_per_seg); 
end

end