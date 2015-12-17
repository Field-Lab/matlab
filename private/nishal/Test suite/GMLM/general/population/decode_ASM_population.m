function [Xdecode,errMap] = decode_ASM_population(fitASM,Y,mask,ttf,rho,lambda,X_ref,B_use)


addpath(genpath('/Volumes/Lab/Users/bhaishahster/cvx'))
cvx_setup
%% parameters
K = fitASM.K;
B = fitASM.B;
Ns = size(K,2);
d =size(K,1);
Nc = size(Y,1);
T = size(Y,2);
dt=1/120;

Tker =gpuArray(zeros(T,T));
for itime=1:T
for jtime = max(itime-29,1):itime
    Tker(itime,jtime) = ttf((itime-jtime)+1)/norm(ttf);
end

% % make it circulant matrix
% if(itime<30)
% Tker(itime,end-30+itime:end) = ttf(end:-1:end-30+itime)/norm(ttf);
% end
end
TTker = Tker'*Tker;
%% initialize
u_k = zeros(d,T);
Z_k = zeros(d,T); % prior mean of stimulus statistics .. 
X_k = zeros(d,T);

scale_cells = sum(Y,2);

ttfspec = sum(abs(fft([ttf;zeros(T-length(ttf),1)])).^2);
stim_spec=0;

 idx=0:T-1 ;
figure;
for idim=1:d
stim_spec = stim_spec+ sum(abs(fft(X_ref(idim,:))).^2);
plot(2*idx/T,abs(fft(X_ref(idim,:))));
hold on;
end

stim_spec = stim_spec/d;
scale_spec = stim_spec/ttfspec;

hold on;
plot(2*idx/T,sqrt(scale_spec)*abs(fft([ttf;zeros(T-length(ttf),1)])),'LineWidth',3);
xlabel('freq (* pi)');
%% iterate
togo = true;
epsi = 10^-5;
K=gather(K);B=gather(B);
while togo
% minimize likelihood
[X_kp1,f_total] = minL_Decode(K,B,Y,rho,Z_k,u_k,X_k,dt,scale_cells);

% project on feasible set
% close all;
% Z_kp1 = projectFeasSet_tf_decode(ttf,X_kp1,u_k,rho);
% Z_kp1 = project_powerSpect_tf_decode(ttf,X_kp1,u_k,rho,lambda,scale_spec);


% another way to incorporate kernels.
tic;
Z_kp1 = gpuArray(zeros(d,T)); u_k = gpuArray(u_k);X_kp1 = gpuArray(X_kp1);
for idim=1:d
Z_kp1(idim,:) = (rho*Tker*((lambda*eye(T,T) + rho*TTker)\(Tker'*(X_kp1(idim,:)+u_k(idim,:))')))';
end
toc;
% update u_k
u_k = u_k + X_kp1 - Z_kp1;

X_k = X_kp1;
Z_k = Z_kp1;

u_k = gather(u_k);
Z_k = gather(Z_k);
X_k = gather(X_k);
[f_total,norm(X_k(:) - Z_k(:))]

if(norm(X_k(:) - Z_k(:)) <epsi)
    togo=false;
    Xdecode= Z_kp1;
end

end

Xdecode=gather(Xdecode);
%%  Compare with reference
Xerr = Xdecode - X_ref;
errMap =zeros(d,1);

for idim=1:d
errMap(idim)  = R_2_value(Xdecode(idim,:)',X_ref(idim,:)');
end

mm = reshape(1:1600,[40,40]);
iidx = mm(mask);

[r,c] = find(repelem(reshape(mask,[40,40])',20,20));
  u=mean(errMap)*ones(1600,1);
    u(iidx) = errMap;
    xx = reshape(u,[40,40])';
    xx = repelem(xx,20,20);
    
    figure;
    imagesc(xx(min(r):max(r),min(c):max(c)));
    axis image
    set(gca,'xTick',[]);
    set(gca,'yTick',[]);
    colormap gray
   
    errImage=xx;
    hold on;
    cols = 'rgbm';
    ncell = (size(Y,1));
    for icell=1:ncell
    plot(B_use{icell}(:,1)-min(c)+1,B_use{icell}(:,2)-min(r)+1,'Color',cols(icell),'LineWidth',3);
    end
end