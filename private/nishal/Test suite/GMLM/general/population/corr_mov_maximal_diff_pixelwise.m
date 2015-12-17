function [corr_pixa,corr_pixb] = corr_mov_maximal_diff_pixelwise(mov_a,mov_b,mov,pc)

corr_pixa=zeros(size(mov_a,1),size(mov_a,2));
corr_pixb=zeros(size(mov_b,1),size(mov_b,2));

for idim1=1:size(mov_a,1)
    for idim2 = 1:size(mov_a,2)
        seqa = squeeze(mov_a(idim1,idim2,:));
        seqb = squeeze(mov_b(idim1,idim2,:));
        diff = abs(seqa-seqb);
        thr = prctile(diff,pc);
        idx = diff>thr;
        seq = squeeze(mov(idim1,idim2,:));
         
        corr_pixa(idim1,idim2) = corr(seqa(idx),seq(idx));
        corr_pixb(idim1,idim2) = corr(seqb(idx),seq(idx));
    end
end

end