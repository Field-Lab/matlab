function STA = calc_sta_mov(movie,spks,usrDepth)
dims = size(movie);
STA =zeros(dims(1),dims(2),usrDepth);

idx = 1:length(spks);
iidx = idx(spks>0 & idx>usrDepth);

for itime = iidx
STA = STA  + movie(:,:,itime:-1:itime-usrDepth+1)*spks(itime);
end
STA = STA/length(iidx);


end