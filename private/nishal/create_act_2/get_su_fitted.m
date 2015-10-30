function [stas_t,stas_r] = get_su_fitted(folder,cellIDs,nSUs)

stas_t = cell(1,1);
stas_r = cell(1,1);

jfilt=1;icell=1;
for cellID=cellIDs
   data = load(['/Volumes/Lab/Users/bhaishahster/',folder,sprintf('/Cell_%d.mat',cellID)]);
     nSU=nSUs(icell);
    st = data.fitGMLM_log{nSU}.stas_true;
    sr = data.fitGMLM_log{nSU}.stas_true;
   
    for ifilt = 1:nSU
    if(data.fitGMLM_log{nSU}.good_su_idx(ifilt)==1)
        stas_t{jfilt} = st{ifilt};
        stas_r{jfilt} = sr{ifilt};
        jfilt=jfilt+1;
    end
    end
    
    icell=icell+1;
end

end