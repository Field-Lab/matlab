function [mov_fit2,binnedSpikeResponses_tr_use]=cleanup_segments_btwn_saccades(mov_fit,binnedSpikeResponses,cut_time,nSeg,cut_edge)

fr_per_seg = size(mov_fit,2)/nSeg;
new_fr_per_seg = (fr_per_seg - cut_time);
mov_fit2 = zeros(size(mov_fit,1),new_fr_per_seg);

remove_seconds = ceil(cut_time/120);

icnt=1;

for imov= 1:nSeg
    for isec= remove_seconds+1:fr_per_seg/120  
    mov_fit2(:,icnt:icnt + 120 - 2*cut_edge) = mov_fit(:,(imov-1)*fr_per_seg + (isec-1)*120 + cut_edge   :(imov-1)*fr_per_seg + isec*120 - cut_edge);
    binnedSpikeResponses_tr_use(icnt:icnt + 120 - 2*cut_edge) = binnedSpikeResponses((imov-1)*fr_per_seg + (isec-1)*120 + cut_edge   :(imov-1)*fr_per_seg + isec*120 - cut_edge); 
    icnt = icnt+ 120-2*cut_edge;
    end
end




end