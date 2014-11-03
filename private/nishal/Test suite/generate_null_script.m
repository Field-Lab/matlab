
if(movie_idx==1)
mov_gen2(319-310,158-150,:)=0.5; % Stimulate a sub-unit
mov_gen2(320-310,160-150,:)=0.5;
mov_gen2=mov_gen2+rand(Filtdim1,Filtdim2,movieLen)*0.3;
end

if(movie_idx==2)
mov_gen2=(double(rand(Filtdim1,Filtdim2,movieLen)>0.5)-0.5)*(0.48/0.50);
end

if(movie_idx==3)
    latency=30;
mov_short=double(rand(Filtdim1,Filtdim2,ceil(movieLen/latency)+2)>0.5)-0.5;
mov_gen2=zeros(Filtdim1,Filtdim2,movieLen);
    for iframe=1:movieLen
        mov_gen2(:,:,iframe)=mov_short(:,:,floor(iframe/latency)+1);
    end
end
if(movie_idx==4)
    latency=20;
mov_short=double(rand(Filtdim1,Filtdim2,ceil(movieLen/latency)+2));
mov_gen2=zeros(Filtdim1,Filtdim2,movieLen);
    for iframe=1:movieLen
        mov_gen2(:,:,iframe)=mov_short(:,:,floor(iframe/latency)+1);
    end
end

if(movie_idx==5)

no_mov=10;%10;
start_image=20;

movie_time=120*no_mov;%size(movie,3);


movie_log=zeros(320,160,movie_time);
icnt=0;
for imov=start_image:no_mov+start_image-1
    imov
    icnt=icnt+1;
movv=load(sprintf('/Volumes/Data/stimuli/movies/eye-movement/NSbrownian_code/matfiles/movie_chunk_%d.mat',imov));
movv=movv.movie;
movie_log(:,:,(icnt-1)*120+1:icnt*120)=movv/255;
end
movie2=movie_log;
clear movie_log
var64=Filtdim1;
var32=Filtdim2;
mov=zeros(var64,Filtdim2,movie_time);
for itime=1:movie_time
mov(:,:,itime)=(imresize(movie2(:,:,itime),[var64,var32],'bilinear','Antialiasing',true)-0.5); % Doubt!!
%mov(:,:,itime)=movie2(var64,32,itime);
display('Resizing a snippet');
end

mov_gen2=mov;
end

mov_gen2(:,:,1:120)=0;


% Compute null space movie
addpath('../create_act_2/');

mov2=zeros(Filtdim1,Filtdim2,movieLen+2*120);
mov2(:,:,121:end-120)=mov_gen2;
stas=cell(1,1);
stas{1}=zeros(Filtdim1,Filtdim2,1,Filtlen);
for itime=1:Filtlen
    itime;
%stas{1}(:,:,1,itime)=STA(:,:,itime); % TODO use B/W STA
    if(sta_null==1)
    stas{1}(:,:,1,itime) = squeeze(STA(:,:,itime));
    end
    
    if(sta_null==2)
    stas{1}(:,:,1,itime) = squeeze(sum(full_fit(:,:,:,itime),3))';
    end
end

 [mov_orig2,mov_new2]=fourier_project_2(stas,mov2);
 
 % See movies
 figure;
 for itime=400:410%movieLen
     itime
  subplot(2,1,1);
  imagesc(mov_orig2(:,:,itime));
  caxis([min(mov_orig2(:)),max(mov_orig2(:))]);
  colormap gray
  colorbar
  axis image
 
  subplot(2,1,2);
  imagesc(mov_new2(:,:,itime));
  caxis([min(mov_new2(:)),max(mov_new2(:))]);
  colormap gray
  colorbar
  axis image
 
  pause(1/120)
 end
 
 % Movie histograms + other stuff 
 figure;
 % Movie histograms 
 subplot(2,1,1);
 [C,X]=hist(mov_orig2(:),100);
 [C2,X2]=hist(mov_new2(:),100);
 plotyy(X,C,X2,C2);
 legend('Original','Null space');
 title('Movie pixel histograms');
 
 % Discretize the movie 
 mov_params.mov_type='bw';
 mov_params.deviation=0.48;
 mov_params.scaling_loss=0.03;

[mov_orig,mov_modify_new]=movie_post_process(mov_orig2,mov_new2,mov_params);
 % Movie histograms 
 subplot(2,1,2);
 [C,X]=hist(mov_orig(:),100);
 [C2,X2]=hist(mov_modify_new(:),100);
 plotyy(X,C,X2,C2);
 legend('Original','Null space');
 title(' Post-processed Movie pixel histograms');
 
 % Discretize ???? 
dum_mov=mov_modify_new;
dum_mov = round((dum_mov)*(127.5/0.5)+127.5);
dum_mov(dum_mov>255)=255;
dum_mov(dum_mov<0)=0;
dum_mov=(dum_mov-127.5)*(0.5/127.5);
mov_modify_new=dum_mov;

dum_mov=mov_orig;
dum_mov = round((dum_mov)*(127.5/0.5)+127.5);
dum_mov(dum_mov>255)=255;
dum_mov(dum_mov<0)=0;
dum_mov=(dum_mov-127.5)*(0.5/127.5);
mov_orig=dum_mov;