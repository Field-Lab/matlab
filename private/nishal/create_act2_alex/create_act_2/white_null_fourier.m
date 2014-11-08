% Try white noise generation in null space using Fourier transform methods
load('/Volumes/Analysis/nishal/Demo/null_data.mat');
ssta=squeeze(stas{3});
ssta=ssta(end:-1:1,end:-1:1,:);


movie_len=120;
sta_mov=zeros(64,32,movie_len);
sta_mov(:,:,1:30)=ssta;
tic;
sta_mov_ft=fftn(sta_mov);
toc;


%mov_ft_re= randn(size(sta_mov_ft)) ; % Sample white noise movie in Fourier domain ? 
%mov_ft_im = randn(size(sta_mov_ft));
%mov_ft = mov_ft_re + sqrt(-1)*mov_ft_im;
%%mov_orig=0*sta_mov;
%%mov_orig(:,:,31:60) = squeeze(stas{1});%randn(size(sta_mov_ft))*120;
mov_orig = randn(size(sta_mov_ft))*120;
%mov_orig(:,:,1:30)=0;
mov_orig(:,:,end-50:end)=0;
mov_ft=fftn(mov_orig);  % In the end, we want to reduce the cost of taking this FFT by sampling a real white movie in frequency space!
%mov_ft_conj= conj(mov_ft(end:-1:1,end:-1:1,end:-1:1));
mov_ft_conj = conjugate_transpose_mov(mov_ft);

%mov_ft =0.5*( mov_ft + conj(mov_ft(end:-1:1,end:-1:1,end:-1:1)));
% Real part is symmetric, Conjugate part is antisymmetric
mov_ft_re = real(mov_ft);
mov_ft_im = imag(mov_ft);
%mov_orig = ifftn(mov_ft);
%TODO - FFT v/s DFT stuff

% Filter movie 
mov_new_ft_re=zeros(size(sta_mov_ft));
mov_new_ft_im=zeros(size(sta_mov_ft));

for iftime=movie_len:-1:1
    iftime
    s_re = real(sta_mov_ft(:,:,iftime)); s_re=s_re(:);
    s_im = imag(sta_mov_ft(:,:,iftime)); s_im=s_im(:);
    A=[s_re',-s_im';
        s_im',s_re'];
    
    m_real = mov_ft_re(:,:,iftime);
    m_im = mov_ft_im(:,:,iftime);
    m=[m_real(:);m_im(:)];
    
   % null_A = (eye(size(A,2)) - A'*((A*A')\A));
    m=m - A'*((A*A')\(A*m));
    m_new_real = reshape(m(1:length(m)/2),size(m_real));
    m_new_im = reshape(m(length(m)/2+1:end),size(m_im));
    
    mov_new_ft_re(:,:,iftime) = m_new_real;
    mov_new_ft_im(:,:,iftime) = m_new_im;
    norm(A*m)
end

% Fill remaining movie time by conjugate symmetry ..
mov_new_ft = mov_new_ft_re+sqrt(-1)*mov_new_ft_im;
mov_new_ft_conj = conjugate_transpose_mov(mov_new_ft);
mov_new_ft_r = 0.5*(mov_new_ft+mov_new_ft_conj);
%

%
mov_new=ifftn(mov_new_ft_r,'symmetric');

%% Verify in fourier domain

lx=[];
for iftime=1:movie_len
    iftime
    s_re = real(sta_mov_ft(:,:,iftime)); s_re=s_re(:);
    s_im = imag(sta_mov_ft(:,:,iftime)); s_im=s_im(:);
    A=[s_re',-s_im';
        s_im',s_re'];
    
    m_real =real(mov_new_ft_r(:,:,iftime));%real(mov_new_ft_r(:,:,iftime));%mov_ft_re(:,:,iftime);
    m_im = imag(mov_new_ft_r(:,:,iftime));%imag(mov_new_ft_r(:,:,iftime));% mov_ft_im(:,:,iftime);
    m=[m_real(:);m_im(:)];
    
   % null_A = (eye(size(A,2)) - A'*((A*A')\A));
  
    op=A*m;
    lx=[lx;op(1)+sqrt(-1)*op(2)];
%lx=[lx;norm(A*m)];
end
plot((ifft(lx,'symmetric')))
%%

figure;
for itime=1:movie_len
    subplot(2,1,1);
    imagesc((mov_orig(:,:,itime)));
colormap gray
caxis([min((mov_orig(:))),max((mov_orig(:)))]);
colorbar
%title(sprintf('Time: %d',itime));
axis image

    subplot(2,1,2);
imagesc((mov_new(:,:,itime)));
colormap gray
caxis([min((mov_orig(:))),max((mov_orig(:)))]);
colorbar
axis image
%title(sprintf('Time: %d',itime));
pause(1/120);
end

%% Cell Response

figure;

%cell_resp = Ax(xs,mov_orig,1200,1);
sz=max(size(mov_orig,3)-size(stas{1},4) + 1, 0);
cell_resp =reshape(convn(mov_orig,ssta,'valid'),[sz,1]);
plot(cell_resp,'r');hold on
sz=max(size(mov_orig,3)-size(stas{1},4) + 1, 0);
cell_resp =reshape(convn(mov_new,ssta,'valid'),[sz,1]);
plot(cell_resp);

%%
figure;
for itime=1:movie_len
subplot(3,1,1);
imagesc(imag(mov_ft(:,:,itime)));
colorbar
subplot(3,1,2)
imagesc(imag(mov_ft_conj(:,:,itime)));
colorbar
subplot(3,1,3)
imagesc(imag(mov_ft(:,:,itime)-mov_ft_conj(:,:,itime)));
colorbar

pause(0.5)
end

%%
v=convn(sta_mov,mov_new);


plot(squeeze(v(64,32,120-91+1:120)));
hold on
sz=max(size(mov_orig,3)-size(stas{1},4) + 1, 0);
cell_resp =reshape(convn(mov_new,ssta,'valid'),[sz,1]);
plot(cell_resp,'r');


ft1 = fftn(sta_mov);
ft2=fftn(mov_new);
x=ft1.*ft2;
l = zeros(size(x,3),1);
for i=1:length(l)
l(i)=sum(sum(x(:,:,i)));
end
plot(abs(l))

plot(abs(ifftn(l)))