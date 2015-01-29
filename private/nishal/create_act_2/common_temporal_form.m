function stas_sp_current = common_temporal_form(stas,CellMasks)

temporal_mean=squeeze(mean(mean(mean(stas,3),1),2));
       [r,c]=find(CellMasks==1);
        stas_sp_current=zeros(size(stas,1),size(stas,2));
       for ipix=1:length(r)
         stas_sp_current(r(ipix),c(ipix))=(squeeze(sum(stas(r(ipix),c(ipix),1,:),3))'*temporal_mean)/(temporal_mean'*temporal_mean);
       end
end