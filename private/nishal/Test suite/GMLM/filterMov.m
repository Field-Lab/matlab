function maskedMov= filterMov(mov,mask,time_c)
%% 
Filtlen=length(time_c);
Filtdim1=size(mov,1);
Filtdim2=size(mov,2);
movieLen=size(mov,3);

mov2=zeros(Filtdim1 ,Filtdim2,movieLen+Filtlen-1);
mov2(:,:,Filtlen:movieLen+Filtlen-1)=mov; % Append zeros before the movie

mov3=zeros(6,6,3,size(mov,3)+Filtlen-1);
mov3(:,:,1,:)=mov2;
mov3(:,:,2,:)=mov2;
mov3(:,:,3,:)=mov2;

mov=mov3;
%%
CellMasks{1}=mask;    

%% 

%time_c(15:end)=0;

figure;
plot(time_c);

mov=squeeze(mean(mov,3));
Filtlen=length(time_c);
filtMov=zeros(size(mov,1),size(mov,2),size(mov,3)-Filtlen+1);
    for x=1:size(mov,1)
        for y=1:size(mov,2)
        filtMov(x,y,:)=conv(squeeze(mov(x,y,:)),time_c,'valid');
        end
    end


%% Filter movie using mask 

maskedMov=zeros(sum(CellMasks{1}(:)),size(filtMov,3));

    for itime=1:size(filtMov,3)
     xx=filtMov(:,:,itime);
     maskedMov(:,itime) = xx(logical(CellMasks{1}));
    end

end