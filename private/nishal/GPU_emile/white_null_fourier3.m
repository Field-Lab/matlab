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

ssta=cell(n_cell,1);
sta_mov=cell(n_cell,1);
sta_mov_ft=cell(n_cell,1);

for icell=1:n_cell
    icell
ssta{icell}=squeeze(stas{icell});
%ssta=zeros(64,32,30);
%ssta(20:40,10:20,1)=1;
ssta{icell}=ssta{icell}(end:-1:1,end:-1:1,:);
sta_mov{icell}=zeros(64*2-1,32*2-1,movie_len);
sta_mov{icell}(1:64,1:32,1:30)=ssta{icell};
tic;
sta_mov_ft{icell}=fftn(sta_mov{icell}).*shift_mat;
toc;
end

%mov_ft_re= randn(size(sta_mov_ft)) ; % Sample white noise movie in Fourier domain ? 
%mov_ft_im = randn(size(sta_mov_ft));
%mov_ft = mov_ft_re + sqrt(-1)*mov_ft_im;
%%mov_orig=0*sta_mov;
%%mov_orig(:,:,31:60) = squeeze(stas{1});%randn(size(sta_mov_ft))*120;
mov_orig = zeros(size(sta_mov_ft{1}));
%mov_orig (64:end,32:end,40:69)= squeeze(stas{3}(:,:,1,:));%randn(64,32,movie_len);
%mov_orig(80:90,40:50,80:160)=10;
mov_orig (64:end,32:end,:)= randn(64,32,movie_len)*120;
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
mov_new_ft_re=zeros(size(sta_mov_ft{1}));
mov_new_ft_im=zeros(size(sta_mov_ft{1}));

tic;
A=zeros(2*n_cell,numel(sta_mov_ft{1}(:,:,1))*2);
for iftime=movie_len:-1:1
   %iftime
    
   
   s_re=cell(2*n_cell,1);
   s_im=cell(2*n_cell,1);
    for icell=1:n_cell*2
   
        if(mod(icell,2)==1)
            s_re{icell} = real(sta_mov_ft{floor((icell-1)/2)+1}(:,:,iftime)); s_re{icell}=s_re{icell}(:);
            s_im{icell} = imag(sta_mov_ft{floor((icell-1)/2)+1}(:,:,iftime)); s_im{icell}=s_im{icell}(:);
            A(icell,:) =[s_re{icell}',-s_im{icell}'];
        else
            A(icell,:)=[s_im{icell-1}',s_re{icell-1}'];
        end
    end
    m_real = mov_ft_re(:,:,iftime);
    m_im = mov_ft_im(:,:,iftime);
    m=[m_real(:);m_im(:)];
    
   % null_A = (eye(size(A,2)) - A'*((A*A')\A));
    m=m - A'*((A*A')\(A*m));
    
    mov_new_ft_re(:,:,iftime) = reshape(m(1:length(m)/2),size(m_real));
    mov_new_ft_im(:,:,iftime) = reshape(m(length(m)/2+1:end),size(m_im));
   % norm(A*m)
    
end
toc;

% Fill remaining movie time by conjugate symmetry ..
mov_new_ft = mov_new_ft_re+sqrt(-1)*mov_new_ft_im;
mov_new_ft_conj = conjugate_transpose_mov(mov_new_ft);
mov_new_ft_r = 0.5*(mov_new_ft+mov_new_ft_conj);
%

%
mov_new=ifftn(mov_new_ft_r,'symmetric');
mov_new2 = mov_new(1:64,1:32,:);


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
