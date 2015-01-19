
% Use code as used in master_run ! 

addpath(genpath('../create_act_2/'));
nFilters=numFilters%4;
stas_coll=cell(nFilters,1);
stas_coll{1}=reSTA;
for iFilter=2:nFilters
stas_coll{iFilter}=WN_uSq{iFilter-1};
end

mov_gen2=zeros(Filtdim1,Filtdim2,movieLen);
movie_idx=2;


if(movie_idx==1)
mov_gen2(319-310,158-150,:)=0.5; % Stimulate a sub-unit
mov_gen2(320-310,160-150,:)=0.5;
mov_gen2=mov_gen2+rand(Filtdim1,Filtdim2,movieLen)*0.3;
end

if(movie_idx==2)
mov_gen2=double(rand(Filtdim1,Filtdim2,movieLen)>0.5)-0.5;
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

mov_gen2(:,:,1:30)=0;

figure;
for itime=40:50
    itime
imagesc(mov_gen2(:,:,itime));
colormap gray
colorbar
axis image
caxis([-0.5,0.5]);
pause(0.01)

end


addpath(genpath('../create_act_2/'));

mov2=zeros(Filtdim1,Filtdim2,movieLen+2*120);
mov2(:,:,121:end-120)=mov_gen2;

stas=cell(1,1);
stas{1}=zeros(Filtdim1,Filtdim2,1,Filtlen);
for itime=1:Filtlen
    itime
    for iFilter=1:nFilters
stas{iFilter}(:,:,1,itime)=stas_coll{iFilter}(:,:,itime);%STA(:,:,itime).*mask;% STA from WN run and masked by significant pixels. Assume I have correct mask from significant_stixles calculation%mask*tf(1,1,1,itime);%STA(:,:,itime);
    end
end

% Clip STA
for iFilter=1:nFilters
stas{iFilter}(:,:,:,15:end)=0;
end

[mov_orig2,mov_new2]=fourier_project_2(stas,mov2);
 
 % See movies
 figure;
 for itime=1000:1000
  subplot(3,2,1);
  imagesc(mov_orig2(:,:,itime));
  caxis([min(mov_new2(:)),max(mov_new2(:))]);
  colormap gray
  colorbar
  axis image
  title('Unmasked original movie');
  
  subplot(3,2,2);
  imagesc(mov_new2(:,:,itime));
  caxis([min(mov_new2(:)),max(mov_new2(:))]);
  colormap gray
  colorbar
  axis image
  title('Unmasked null movie');
  
  subplot(3,2,3);
  imagesc(mov_orig2(:,:,itime).*mask);
  caxis([min(mov_new2(:)),max(mov_new2(:))]);
  colormap gray
  colorbar
  axis image
  title('Masked original movie');  
  
  subplot(3,2,4);
  imagesc(mov_new2(:,:,itime).*mask);
  caxis([min(mov_new2(:)),max(mov_new2(:))]);
  colormap gray
  colorbar
  axis image
  title('Masked null movie');
  
  subplot(3,2,5);
  imagesc((mov_new2(:,:,itime)-mov_orig2(:,:,itime)));
  colormap gray
  colorbar
  axis image
  title('Difference between original and null movies (masked)');
  
  pause(1/120)
 end
 
 % Movie histograms 
  
  mov_gen2_masked=mov_gen2.*repmat(mask,[1,1,size(mov_gen2,3)]);
  mov_new2_center=mov_new2(:,:,121:end-120);
  mov_new2_masked = mov_new2_center.*repmat(mask,[1,1,size(mov_new2_center,3)]);
  
 [C,X]=hist(mov_gen2_masked(mov_gen2_masked~=0),100);
 [C2,X2]=hist(mov_new2_masked(mov_gen2_masked~=0),100);
 figure;
 subplot(2,1,1);
 [ax,h1,h2]=plotyy(X,C,X2,C2);
 xlim(ax(1),[min(X2),max(X2)]);
 xlim(ax(2),[min(X2),max(X2)]);
 
 legend('Original','Null space');
 title('Movie pixel histograms in mask');
 
subplot(2,1,2);

difference =  mov_gen2_masked - mov_new2_masked;
hist(difference(mov_gen2_masked~=0),20);
title('Pixel value difference in mask');

% Subunit Pixel time course
[r,c]=find(logical(mask));
pixCourse=[];
for iPix=1:length(r)
pixCourse(iPix,:)=mov_new2_masked(r(iPix),c(iPix),:);
end

figure;
subplot(2,1,1);
plot(pixCourse(:,500:600)');
title('Subunit Pixel value time courses');

subplot(2,1,2);
plot(max((pixCourse(:,500:700)))','b');
hold on;
plot(min((pixCourse(:,500:700)))','r');
hold on
plot(max(abs(pixCourse(:,500:700)))','k');
hold on
plot([1,200],[0.5,0.5],'m');
hold on;
plot([1,200],[-0.5,-0.5],'m');
legend('max','min','Absolute max','Original movie')
title('Max and Min pixel value accross frames')

