
function movout = spatial_filter(movin,stas,spatial_f_cutoff_lower,spatial_f_cutoff_upper,s1,s2)

%% spatial filtering
%spatial_f_cutoff = 0.5;
% s1=32;s2=16; ??
T  = size(movin,3);

% find frequencies relevant for STA
fft_xx = zeros(s1,s2);
%stas =stas_celltype1;
for icell=1:length(stas)
xx = reshape(stas{icell},[s1,s2]);

% [r,c]=find(xx==max(xx(:)));
% 
% zz=zeros(s1,s2);
% zz
fft_xx = fft_xx + abs(fft2(xx));
end


fft_sp_mask = double( fft_xx<=spatial_f_cutoff_upper*max(fft_xx(:)) &fft_xx>=spatial_f_cutoff_lower*max(fft_xx(:))); % fft_xx<0.5*max(fft_xx(:)) &
figure('Color','w');
subplot(2,1,1);
imagesc(fft_xx');
axis image
title('Average Fourier tranform(magnitude) of RFs')
subplot(2,1,2);
imagesc(fft_sp_mask')
axis image 
title('Selected frequencies')

movout =zeros(s1,s2,T);
 for itime = 1:T
 %xx =  reshape(movin(:,itime),[s1,s2]);
 xx = movin(:,:,itime);
 movout(:,:,itime) = ifft2(fft2(xx).*fft_sp_mask);
 end
 movout = movout * norm(movin(:)) / norm(movout(:));

end