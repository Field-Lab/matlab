function maskedMov= filterMov(mov,mask,time_c)
%% 
Filtlen=length(time_c);
Filtdim1=size(mov,1);
Filtdim2=size(mov,2);
movieLen=size(mov,3);

mov2=zeros(Filtdim1 ,Filtdim2,movieLen+Filtlen-1);
mov2(:,:,Filtlen:movieLen+Filtlen-1)=mov; % Append zeros before the movie

%%
CellMasks{1}=mask;    

 %% Filter movie using mask 
filtMov=mov2;
maskedMov=zeros(sum(CellMasks{1}(:)),size(filtMov,3));

    for itime=1:size(filtMov,3)
     xx=filtMov(:,:,itime);
     maskedMov(:,itime) = xx(logical(CellMasks{1}));
    end

%%
%time_c(15:end)=0;

figure;
plot(time_c);

mov=maskedMov;
Filtlen=length(time_c);
filtMov=zeros(size(mov,1),size(mov,2)-Filtlen+1);
    for x=1:size(mov,1)
        filtMov(x,:)=conv(squeeze(mov(x,:)),time_c,'valid');
    end

maskedMov=filtMov;

end