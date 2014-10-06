% Try white noise generation in null space using Fourier transform methods
load('null_data.mat');


movie_len=1200;
shift_mat=zeros(2*64-1,2*32-1,movie_len);
for ix=1:64*2-1
    for iy=1:32*2-1
        for itime=1:movie_len
            shift_mat(ix,iy,itime) = exp(sqrt(-1)*( (2*pi*(ix-1)*63/(2*64-1)) +(2*pi*(iy-1)*31/(2*32-1)) ) ); 
        end
    end
end
%shift_mat=gpuArray(shift_mat);

ssta=zeros(64,32,30,n_cell);
sta_mov=zeros(64*2-1,32*2-1,movie_len,n_cell);
sta_mov_ft=zeros(size(sta_mov));

%stas=gpuArray(stas);
tic;
for icell=1:n_cell
ssta(:,:,:,icell)=(squeeze(stas{icell}));
%ssta=zeros(64,32,30);
%ssta(20:40,10:20,1)=1;
ssta(:,:,:,icell)=ssta(end:-1:1,end:-1:1,:,icell);
%sta_mov{icell}=zeros(64*2-1,32*2-1,movie_len);
sta_mov(1:64,1:32,1:30,icell)=ssta(:,:,:,icell);
sta_mov_ft(:,:,:,icell)=fftn(squeeze(sta_mov(:,:,:,icell))).*shift_mat;
end
toc;

%mov_ft_re= randn(size(sta_mov_ft)) ; % Sample white noise movie in Fourier domain ? 
%mov_ft_im = randn(size(sta_mov_ft));
%mov_ft = mov_ft_re + sqrt(-1)*mov_ft_im;
%%mov_orig=0*sta_mov;
%%mov_orig(:,:,31:60) = squeeze(stas{1});%randn(size(sta_mov_ft))*120;
mov_orig = zeros(64*2-1,32*2-1,movie_len);
%mov_orig (64:end,32:end,40:69)= squeeze(stas{3}(:,:,1,:));%randn(64,32,movie_len);
%mov_orig(80:90,40:50,80:160)=10;
mov_orig (64:end,32:end,:)= double(rand(64,32,movie_len)>0.5);
mov_orig(1:63,1:31,:)=0; % make sure its 0 where STAs there ?

mov_orig(:,:,1:30)=0;
mov_orig(:,:,end-50:end)=0;

mov_ft=fftn(mov_orig).*shift_mat;  % In the end, we want to reduce the cost of taking this FFT by sampling a real white movie in frequency space!
mov_orig2 = mov_orig(64:end,32:end,:);

%mov_ft_conj= conj(mov_ft(end:-1:1,end:-1:1,end:-1:1));
%mov_ft_conj = conjugate_transpose_mov(mov_ft);

%mov_ft =0.5*( mov_ft + conj(mov_ft(end:-1:1,end:-1:1,end:-1:1)));
% Real part is symmetric, Conjugate part is antisymmetric
mov_ft_re = real(mov_ft);
mov_ft_im = imag(mov_ft);
%mov_orig = ifftn(mov_ft);
%TODO - FFT v/s DFT stuff

% Filter movie 
mov_new_ft_re= zeros(64*2-1,32*2-1,movie_len);
mov_new_ft_im= zeros(64*2-1,32*2-1,movie_len);

tic;

parfor iftime=1:movie_len
   %iftime
  [mov_new_ft_re(:,:,iftime), mov_new_ft_im(:,:,iftime)]= proj_null_fourier(sta_mov_ft(:,:,iftime,:),n_cell,mov_ft_re(:,:,iftime),mov_ft_im(:,:,iftime));
    
end
toc;

% Fill remaining movie time by conjugate symmetry ..
mov_new_ft = mov_new_ft_re+sqrt(-1)*mov_new_ft_im;
mov_new_ft=gather(mov_new_ft);
mov_new_ft_conj = conjugate_transpose_mov(mov_new_ft);
mov_new_ft_r = 0.5*(mov_new_ft+mov_new_ft_conj);
%

%
mov_new_ft_r=(mov_new_ft_r);
mov_new=ifftn(mov_new_ft_r,'symmetric');
mov_new2 = mov_new(1:64,1:32,:);

mov_new2=gather(mov_new2);
mov_orig2=gather(mov_orig2);

%%

cell_resp_o = Ax(stas,mov_orig2,movie_len,n_cell);
cell_resp_m = Ax(stas,mov_new2,movie_len,n_cell);

figure;
subplot(2,1,1);
plot(cell_resp_o);
subplot(2,1,2);
plot(cell_resp_m);

%% 
figure;
imagesc(mean(mov_new2,3));
caxis([min(mov_new2(:)),max(mov_new2(:))]);

colormap gray
colorbar
