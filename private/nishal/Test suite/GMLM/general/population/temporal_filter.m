
function movout = temporal_filter(movin,ttf, temp_f_cutoff)
%%  temporal frequency filtering

%temp_f_cutoff=0.7;

T = size(movin,2);
figure('Color','w');
ftt_total = abs(fft(ttf,T));
plot([0:2/T:2*(T-1)/T],abs(ftt_total))
hold on;

fft_mask = double(ftt_total > temp_f_cutoff*max(ftt_total));
plot([0:2/T:2*(T-1)/T],fft_mask);
legend('ttf',sprintf('Selected mask (%0.02d>max)',temp_f_cutoff));


% filter reconstructed movie
movout = 0*movin;
for idim = 1:size(movin,1)
   movout(idim,:)= ifft(fft(movin(idim,:),T).*fft_mask');
end

movout = movout * norm(movin(:)) / norm(movout(:));

end

