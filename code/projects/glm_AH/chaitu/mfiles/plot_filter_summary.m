function plot_filter_summary(pstar,basepars,neff,dt)


fs = 12;

[K PS CP KSQ] = get_model_filters(pstar,1,basepars,neff);

subplot(2,3,1)
imagesc(K), axis image, set(gca,'FontSize',fs), title(sprintf('Stim filt: b=%0.2f knorm=%0.2f',pstar(get_pars_idx(basepars,1,neff,'b')),norm(K,'fro').*sqrt(dt))), colorbar;

[s1 t1] = get_spacetime_approx(K,1,0);

1;

subplot(2,3,2)
imagesc(reshape(s1,sqrt(basepars.n),sqrt(basepars.n))), axis image, colorbar;
set(gca,'FontSize',fs), title('Rank 1 spatial profile')

subplot(2,3,3)
cla, plot(makeaxis(dt,basepars.Mk),t1), set(gca,'FontSize',fs), title('Rank 1 temporal profile');
hold on, plot(makeaxis(dt,basepars.Mk),zeros(basepars.Mk,1),'k--','LineWidth',2);

subplot(2,3,4)
cla, plot(makeaxis(dt/basepars.fac,basepars.Mhist),PS), set(gca,'FontSize',fs), title('Postspike filter: long timescale'), xlabel('time (s)'), ylim([-0.1 0.1]);
hold on, plot(makeaxis(dt/basepars.fac,basepars.Mhist),zeros(basepars.Mhist,1),'k--','LineWidth',2);

subplot(2,3,5)
cla, plot(makeaxis(dt/basepars.fac,basepars.Mhist),PS), set(gca,'FontSize',fs), title('Postspike filter: short timescale'), xlabel('time (s)'), ylim([-10 1]), xlim([0 0.1]);
hold on, plot(makeaxis(dt/basepars.fac,basepars.Mhist),zeros(basepars.Mhist,1),'k--','LineWidth',2);

if (neff > 1 && ~isempty(CP) && nargin > 5)
    subplot(2,3,6)
    hold on , plot (makeaxis(dt/basepars.fac,basepars.Mcoup),CP), set(gca,'FontSize',fs), title('Coupling filters');
    hold on, plot(makeaxis(dt/basepars.fac,basepars.Mhist),zeros(size(psaxis)),'k--','LineWidth',2);
end

if (~isempty(KSQ))
    subplot(2,3,6)
    imagesc(KSQ), axis image, set(gca,'FontSize',fs), title(sprintf('Stim sq filt: knorm=%0.2f',norm(KSQ,'fro').*sqrt(dt))), colorbar;
end

if (isfield(basepars,'ext_timepts') && ~isempty(basepars.ext_timepts))
    subplot(2,3,6);
    plot_ext_signal(pstar,basepars,dt,neff);
    set(gca,'FontSize',fs), title('Extrinsic signal fit'), xlabel('time (s)');
end