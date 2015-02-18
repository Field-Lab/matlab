function [mov_orig2,mov_new2]=fourier_project256x256(stas,mov)

display('Starting Fourier Computation');
var64=256;%size(stas{1},1);

movie_time=size(mov,3);
filt_len=size(stas{1},4);
movie_len=movie_time+3*filt_len;
n_cell=length(stas);

tic;
shift_mat=zeros(2*var64-1,2*256-1,movie_len);


% Slower code
% for ix=1:256*2-1
%     for iy=1:256*2-1
%         for itime=1:movie_len
%             shift_mat(ix,iy,itime) = exp(sqrt(-1)*( (2*pi*(ix-1)*(256-1)/(2*256-1)) +(2*pi*(iy-1)*31/(2*256-1)) ) ); 
%         end
%     end
% end

% Faster code
ix=1:256*2-1;
iy = 1:256*2-1;
shift_mat1 = repmat(( (2*pi*(ix-1)*(256-1)/(2*256-1)))',[1,length(iy)]);
shift_mat2 = repmat((2*pi*(iy-1)*(256-1)/(2*256-1)),[length(ix),1]);
shift_mat3=shift_mat1+shift_mat2;
shift_mat3=exp(sqrt(-1)*(shift_mat3));
shift_mat=repmat(shift_mat3,[1,1,movie_len]);



%shift_mat=gpuArray(shift_mat);
toc;

ssta=zeros(var64,256,30,n_cell);
sta_mov=zeros(var64*2-1,256*2-1,movie_len,n_cell);
sta_mov_ft=zeros(size(sta_mov));

%stas=gpuArray(stas);
tic; % VERY slow! Make it faster  ?
for icell=1:n_cell
    icell
    %display('CAREFUL, Using only 1 cell!')
ssta(:,:,:,icell)=(squeeze(stas{icell}));
ssta(:,:,:,icell)=ssta(end:-1:1,end:-1:1,:,icell);
sta_mov(1:var64,1:256,1:30,icell)=ssta(:,:,:,icell);
sta_mov_ft(:,:,:,icell)=fftn(squeeze(sta_mov(:,:,:,icell))).*shift_mat;
end
toc;


mov_orig = zeros(var64*2-1,256*2-1,movie_len);
mov_orig (var64:end,256:end,end-movie_time+1:end)= mov;
mov_orig(1:var64-1,1:256-1,:)=0; % make sure its 0 where STAs there ?

mov_orig(:,:,1:30)=0;
mov_orig(:,:,end-50:end)=0;

mov_ft=fftn(mov_orig).*shift_mat;  % In the end, we want to reduce the cost of taking this FFT by sampling a real white movie in frequency space!
mov_orig2 = mov_orig(var64:end,256:end,end-movie_time+1:end);

% Real part is symmetric, Conjugate part is antisymmetric
mov_ft_re = real(mov_ft);
mov_ft_im = imag(mov_ft);
%TODO - FFT v/s DFT stuff

% Filter movie 
mov_new_ft_re= zeros(var64*2-1,256*2-1,movie_len);
mov_new_ft_im= zeros(var64*2-1,256*2-1,movie_len);

tic;

parfor iftime=1:movie_len
   iftime
  [mov_new_ft_re(:,:,iftime), mov_new_ft_im(:,:,iftime)]= proj_null_fourier256x256(sta_mov_ft(:,:,iftime,:),n_cell,mov_ft_re(:,:,iftime),mov_ft_im(:,:,iftime));
    
end
toc;

% Fill remaining movie time by conjugate symmetry ..
mov_new_ft = mov_new_ft_re+sqrt(-1)*mov_new_ft_im;
mov_new_ft=gather(mov_new_ft);
mov_new_ft_conj = conjugate_transpose_mov(mov_new_ft);
mov_new_ft_r = 0.5*(mov_new_ft+mov_new_ft_conj);
%

%
% Try 
%mov_new_ft_r=mov_new_ft; % DEBUG
%
mov_new=ifftn(mov_new_ft_r); % DEBUG
mov_new2 = mov_new(1:var64,1:256,:);
mov_new2=gather(mov_new2);
mov_new2=mov_new2(:,:,end-movie_time+1:end);

mov_orig2=gather(mov_orig2);

cell_resp_o = Ax(stas,mov_orig2,movie_time,n_cell);
cell_resp_m = Ax(stas,mov_new2,movie_time,n_cell);

figure;
subplot(2,1,1);
plot(cell_resp_o);
subplot(2,1,2);
plot(cell_resp_m);
end
% %%
% 
% cell_resp_o = Ax(stas,mov_orig2,movie_len,n_cell);
% cell_resp_m = Ax(stas,mov_new2,movie_len,n_cell);
% 
% figure;
% subplot(2,1,1);
% plot(cell_resp_o);
% subplot(2,1,2);
% plot(cell_resp_m);
% 
% %% 
% figure;
% imagesc(mean(mov_new2,3));
% caxis([min(mov_new2(:)),max(mov_new2(:))]);
% 
% colormap gray
% colorbar
