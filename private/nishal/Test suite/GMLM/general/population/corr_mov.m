function corr_pix = corr_mov(mov1,mov2)

corr_pix=zeros(size(mov1,1),size(mov1,2));
for idim1=1:size(mov1,1)
    for idim2 = 1:size(mov1,2)
        seq1 = squeeze(mov1(idim1,idim2,:));
        seq2 = squeeze(mov2(idim1,idim2,:));
        corr_pix(idim1,idim2) = corr(seq1,seq2);
    end
end

end