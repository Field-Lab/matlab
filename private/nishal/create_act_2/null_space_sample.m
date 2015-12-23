load('/Volumes/Analysis/nora/NSEM/GLM_Output/rk1_MU_PS_noCP_noSU_p8IDp8/standardparams/WN_mapPRJ/2012-08-09-3/OFFPar_1471.mat');
%%
k_space=fittedGLM.linearfilters.Stimulus.space_rk1;
k_time=fittedGLM.linearfilters.Stimulus.time_rk1;
%k_time=k_time-mean(k_time);
filt_len=length(k_time);
filt_dim=13;
%% 
movie_time=120%200;
mov = rand(filt_dim,filt_dim,movie_time);
%mov=movie;
for itime=1:movie_time
   imagesc(mov(:,:,itime));
   colormap gray
   title(sprintf('Time: %d',itime));
   pause(0.1);
end

%%
k=zeros(filt_dim,filt_dim,movie_time);
for itime=1:length(k_time)
k(:,:,itime)=k_space*k_time(itime,1);
end

%%
filt_mat= sparse(zeros(filt_dim*filt_dim*movie_time));

for itime=1:movie_time
   
    for iitime=1:min(itime,filt_len)
        
    filt_mat(itime,itime*filt_dim*filt_dim-(iitime)*filt_dim*filt_dim+1:itime*filt_dim*filt_dim-(iitime-1)*filt_dim*filt_dim)=k_space(:)'*k_time(iitime);
    
    end
end

%% 
mov_reshape=zeros(filt_dim*filt_dim*movie_time,1);
itime=0;
for i=1:filt_dim*filt_dim:length(mov_reshape)
    itime=itime+1;
    xx=mov(:,:,itime);
  mov_reshape(i:i+filt_dim*filt_dim-1)=  xx(:);
end

%%

dd=filt_dim1*filt_dim1*movie_time;
tic;
cvx_begin
variable mov_modify(dd)
filt_mat*mov_modify==0
minimize (norm(mov_modify-mov_reshape))
cvx_end
toc;
%% 
%filt_proj=(eye(dd)-filt_mat'*inv(filt_mat*filt_mat')*filt_mat);
tic;
mov_modify2=mov_reshape-filt_mat'*(filt_mat*filt_mat')\filt_mat*mov_reshape;
toc;

%%
%%
mov_mod=zeros(filt_dim,filt_dim,movie_time);

for itime=1:movie_time
mov_mod(:,:,itime)=reshape(mov_modify(1+(itime-1)*filt_dim*filt_dim:itime*filt_dim*filt_dim),[filt_dim, filt_dim]);
end

%%
f=figure;
subplot(1,2,1)
imagesc(k_space);colormap(gray)
colorbar
axis square
title('Space filter')

subplot(1,2,2)
plot([-29:0],k_time(end:-1:1));
title('Time filter');

print(f,'Filters_1.png','-dpng');


%%  
[movfile,movpath] = uiputfile('*.avi','Save Movie As');
    movieFileName = [movpath movfile];
    writerObj = VideoWriter(movieFileName);
    writerObj.FrameRate = 3;
    open(writerObj);

f = figure; set(f,'Position',[100 360 1000 550]);
set(f,'Color','white');

for itime=1:movie_time
    
    subplot(2,2,1)
     imagesc(mov(:,:,itime),[0,1]);
   title(sprintf('Time: %d',itime));
   axis square
   
colorbar
subplot(2,2,2)
imagesc(mov_mod(:,:,itime),[0,1]);colormap(gray)
title(sprintf('Time: %d',itime));
axis square
colorbar
 subplot(2,2,3)
imagesc(mov_mod(:,:,itime)-mov(:,:,itime),[-0.63,0.1]);colormap(gray)
title(sprintf('Time: %d',itime));
axis square
colorbar


subplot(2,2,4)
imagesc(k_space);colormap(gray)
colorbar
axis square
title('RF in data')


 M = getframe(f);
 writeVideo(writerObj,M);
 pause(0.5)
end

 close(writerObj);

 % All temporal structure lost!- even in stimulus filter, because of first
 % order structure!
 
 %% 
 % Verify!
 f=figure;
 subplot(1,2,1);
 plot(filt_mat*mov_modify);
 title('Both');
 
 % Make filt_mat2 - see if only in spatial null space? 
 filt_mat2= sparse(zeros(filt_dim*filt_dim*movie_time));

for itime=1:movie_time
   
    for iitime=1:1%min(itime,filt_len)
        
    filt_mat2(itime,itime*filt_dim*filt_dim-(iitime)*filt_dim*filt_dim+1:itime*filt_dim*filt_dim-(iitime-1)*filt_dim*filt_dim)=k_space(:)'*k_time(iitime);
    
    end
end
subplot(1,2,2);
plot(filt_mat2*mov_modify);
title(' only spatial filter');
print(f,'Stimulus_2_inner_prod.png','-dpng');
