clear all
clf
N = 64;
x = -N/2:N/2-1;
sigma_m_f = 10/N;
sigma_q_f = 20/N;

% M - big subunits
clear M
M_sc = normpdf(x,0,1/sigma_m_f);
M_sc = circshift(M_sc(1,:),[0 N/2]);
M_sc = M_sc/sum(M_sc);
for i=2:N;
    M_sc(i,:) = circshift(M_sc(i-1,:),[0 1]);
end

% Q - small subunits
clear Q
Q_sc = normpdf(x,0,1/sigma_q_f);
Q_sc = circshift(Q_sc(1,:),[0 N/2]);
Q_sc = Q_sc/sum(Q_sc);
for i=2:N;
    Q_sc(i,:) = circshift(Q_sc(i-1,:),[0 1]);
end

M_sc = M_sc(1:8:end,:); % 8 or 32
Q_sc = Q_sc(1:4:end,:); % 4 or 16

clrs_m = cmap_hsv_dark(size(M_sc,1));
clrs_q = cmap_hsv_dark(size(Q_sc,1));
subplot(2,1,1);
hold on
for is=1:size(M_sc,1);
    h = plot(M_sc(is,:));
    set(h,'Color',clrs_m(is,:));
    set(h,'LineWidth',3);
end
xlabel('Position');
set(gca,'XTick',[]);
ylabel('Response');
set(gca,'YTick',[]);
xlim([0 N+1]);
ylim([0 max(M_sc(1,:))+0.01]);

subplot(2,1,2);
hold on
for is=1:size(Q_sc,1)
    h = plot(Q_sc(is,:));
    set(h,'Color',clrs_q(is,:));
    set(h,'LineWidth',3);
end
xlabel('Position');
set(gca,'XTick',[]);
ylabel('Response');
set(gca,'YTick',[]);
xlim([0 N+1]);
ylim([0 max(Q_sc(1,:))+0.01]);
set(gcf,'Position',[858 345 560 420]);
%%
T = 60000;
fit.f.type = 'rectLinear';
fit.f.p = [0 1];
fit.g.type = 'rectLinear';
fit.g.p = [0 1];
train_q.X_ct = randn(N,T);
train_m.X_ct = train_q.X_ct;
Y_st_q = evalNonLin(Q_sc*train_q.X_ct,fit.f);
Y_st_m = evalNonLin(M_sc*train_m.X_ct,fit.f);

B_s_q = ones(1,size(Q_sc,1));
B_s_m = ones(1,size(M_sc,1))*2*sqrt(2);

Z_t_q = evalNonLin(B_s_q*Y_st_q,fit.g);
Z_t_m = evalNonLin(B_s_m*Y_st_m,fit.g);
train_q.R_t = poissrnd(Z_t_q);
train_m.R_t = poissrnd(Z_t_m);

mean(train_q.R_t)
mean(train_m.R_t)

%%
fit_m = getSTC(train_m);
fit_q = getSTC(train_q);
%%
clf
set(gcf,'Position',[858 345 400 625]);
%%
clf
subplot(6,2,2);
plot(fit_q.e(1:12),'k.');
xlim([0 13]);
box off
for i=1:5
    subplot(6,2,2+i*2)
    plot(fit_q.V(:,i),'k','LineWidth',2);
    xlim([0 N+1]);
    ylim([-0.5 0.5]);
    box off
    set(gca,'TickDir','out');
    axis off
end

subplot(6,2,1);
plot(fit_m.e(1:12),'k.');
xlim([0 13]);
box off
for i=1:5
    subplot(6,2,1+2*i)
    plot(fit_m.V(:,i),'k','LineWidth',2);
    xlim([0 N+1]);
    ylim([-0.5 0.5]);
    box off
    set(gca,'TickDir','out');
    axis off
end

%%
fft_q = fftshift(abs(fft(Q_sc(1,:))))';
fft_m = fftshift(abs(fft(M_sc(1,:))))';

%%
clf
iCount = 1;
clear r_filt
for i=1:32
        N = 64;
        x = -N/2:N/2-1;
        sigma_1 = i/N;
        foo1 = normpdf(x,0,1/sigma_1);
        foo1 = foo1/sum(foo1);
        ft1 = fftshift(abs(fft(foo1)));
        r_low(:,iCount) = ft1;
        iCount = iCount + 1;
end
%%
T = 1000;
X_ct = randn(N,T);
X_ct_sin = zeros(size(X_ct));
r = r_low(:,1);
for it=1:T
    noise = X_ct(:,it);
    tmpFFT = vector(fftshift(fft2(noise)));
    tmpTheta = angle(tmpFFT);
    tmpAdjFFT = reshape(complex(r.*cos(tmpTheta),r.*sin(tmpTheta)),size(noise,1),size(noise,2));
    out = real(ifft2(ifftshift(tmpAdjFFT),'symmetric'));
    out = out/std(out(:));
    X_ct_sin(:,it) = out;
end
%%
for iFilt=1:size(r_filt,2)
    Y_st_q = evalNonLin(Q_sc*stim(iFilt).X_ct_sin,fit.f);
    Y_st_m = evalNonLin(M_sc*stim(iFilt).X_ct_sin,fit.f);
    Z_t_q = evalNonLin(B_s_q*Y_st_q,fit.g);
    Z_t_m = evalNonLin(B_s_m*Y_st_m,fit.g);
    train_q.sin(iFilt).R_t = poissrnd(Z_t_q);
    train_m.sin(iFilt).R_t = poissrnd(Z_t_m);
end
mean(train_q.sin(1).R_t)
%%
clear mn_q
clear mn_m
clear vr_q
clear vr_m
for iFilt=1:size(r_filt,2)
    mn_q(iFilt) = mean(train_q.sin(iFilt).R_t);
    mn_m(iFilt) = mean(train_m.sin(iFilt).R_t);
    vr_q(iFilt) = var(train_q.sin(iFilt).R_t);
    vr_m(iFilt) = var(train_m.sin(iFilt).R_t);
end