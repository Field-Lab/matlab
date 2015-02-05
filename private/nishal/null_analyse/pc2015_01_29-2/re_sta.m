% Re-STA

%Use spkCondCollformat(icond).spksColl
% Use condMovies
MovLen=size(condMovies{1},1);
re_stas=cell(4,1);


for icond=1:4
    sta_calc=zeros(30,Filtdim1,Filtdim2);
icnt=0;
for itrial=1:30
    for itime=31:MovLen
      if(spkCondCollformat(icond).spksColl(itrial,itime)==1)
      sta_calc=sta_calc+condMovies{icond}(itime-30+1:itime,:,:);
      icnt=icnt+1;
      end
    end
end

re_stas{icond}=sta_calc/icnt;
end


%%
figure;
for itime=30:-1:1

for icond=1:4
subplot(2,2,icond);
imagesc(squeeze(re_stas{icond}(itime,:,:)));
axis image
colormap gray
colorbar
title(sprintf('Time :%d Cond: %d',30-itime+1,icond));
caxis([min(re_stas{icond}(:)),max(re_stas{icond}(:))])
end
pause
end