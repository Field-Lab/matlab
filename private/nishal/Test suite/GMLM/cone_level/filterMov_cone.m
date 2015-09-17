function maskedMov= filterMov_cone(mov,mask,time_c)
%% 
Filtlen=length(time_c);
Filtdim1=size(mov,1);
Filtdim2=1;
movieLen=size(mov,2);

mov2=gpuArray(zeros(Filtdim1 ,movieLen+Filtlen-1));
mov2(:,Filtlen:movieLen+Filtlen-1)=mov; % Append zeros before the movie

mov3=gpuArray(zeros(Filtdim1,1,3,movieLen+Filtlen-1));
mov3(:,1,1,:)=mov2;
mov3(:,1,2,:)=mov2;
mov3(:,1,3,:)=mov2;

mov=mov3;
%%
CellMasks{1}=mask;    

%% 

%time_c(15:end)=0;
% 
% figure;
% plot(time_c);

mov=squeeze(mean(mov,3));
Filtlen=length(time_c);
filtMov=gpuArray(zeros(size(mov,1),size(mov,2)-Filtlen+1));
    for x=1:size(mov,1)
        filtMov(x,:)=conv(squeeze(mov(x,:)),time_c,'valid');
    end


%% Filter movie using mask 

maskedMov=gpuArray(zeros(sum(CellMasks{1}(:)),size(filtMov,2)));

    for itime=1:size(filtMov,2)
     xx=filtMov(:,itime);
     maskedMov(:,itime) = xx(logical(CellMasks{1}));
    end

end