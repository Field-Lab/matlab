function maskedMov = mov_slice_flattened(mov2,mask)

filtMov=mov2;
maskedMov=zeros(sum(mask(:)),size(filtMov,3));

    for itime=1:size(filtMov,3)
     xx=filtMov(:,:,itime);
     maskedMov(:,itime) = xx(logical(mask));
    end

end