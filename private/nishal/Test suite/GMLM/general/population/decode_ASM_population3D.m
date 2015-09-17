function [Xdecode,errMap] = decode_ASM_population3D(fitASM,Y,mask,ttf,rho,lambda,X_ref,B_use)


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
scale_cells = sum(Y,2);

%% iterate
Tker =gpuArray(zeros(T,T));
for jtime=1:T
for itime = jtime:min(jtime+29,T)
    Tker(itime,jtime) = 1000*ttf((itime-jtime)+1);
end
end

optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','iter-detailed',...  %'iter-detailed',... 
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','off');

 Xinit = randn(d,T);
% 
tic;
xfinal = fminunc(@(x)minL_dec3D(x,K,B,Y,dt,d,T,scale_cells,ttf,Tker),Xinit(:),optim_struct);
toc;
Xdecode= xfinal;
Xdecode = reshape(Xdecode,[d,T]);

 Xdecode= filterMov_cone(Xdecode,logical(ones(d,1)),squeeze(ttf));
% 
% x= Xinit;
% for iter=1:2000
%     eta= 0.1/iter;
% [fval,x] = minL_dec3D_gd(x,K,B,Y,dt,d,T,scale_cells,1000*ttf,eta,Tker);
% [fval,norm(x(:))]
% 
% end
% 
% Xdecode = gather(filterMov_cone(x,logical(ones(d,1)),squeeze(ttf)));
      
%%  Compare with reference
Xerr = Xdecode - X_ref;
errMap =zeros(d,1);

for idim=1:d
errMap(idim)  = R_2_value(Xdecode(idim,:)',X_ref(idim,:)');
end

mm = reshape(1:3200,[80,40]);
iidx = mm(mask);

[r,c] = find(repelem(reshape(mask,[80,40])',20,20));
  u=mean(errMap)*ones(3200,1);
    u(iidx) = errMap;
    xx = reshape(u,[80,40])';
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