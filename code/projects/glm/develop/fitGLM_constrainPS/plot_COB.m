% hack plot to show how the pasis change to do constarined searches
% AKHeitman 2015-07-06

close all;
GLMPars           = GLMParams;
bin_size = 1/1200;


basis_params  = GLMPars.spikefilters.ps;
basis_params.filternumber = 10;
ps_basis      = prep_spikefilterbasisGP(basis_params,bin_size);

ps_basis_0 = ps_basis; clear ps_basis
v        = sum(ps_basis_0,1);
v        = v / norm(v) ;
orthog_v = null(v);
COB      = [v', orthog_v] ;
ps_basis = (inv(COB) * ps_basis_0')' ;


figure;
subplot(1,2,1);
plot((ps_basis(:,1)),'k','linewidth', 2); hold on
axis square
legend('only non zero-gain term','location','southeast')
plot((ps_basis),'linewidth',1); hold on
plot((ps_basis(:,1)),'k','linewidth', 2);
subframe_tick = 0:30:120;
msec_tick   = 0:25:100;
set(gca, 'xtick', subframe_tick, 'xticklabel',msec_tick)
xlim([0 120])

xlabel('milliseconds');
title('Orthonormal Change of Basis')

subplot(1,2,2);
plot(ps_basis_0,'linewidth',1); hold on; axis square;
xlabel('milliseconds');
subframe_tick = 0:30:120;
msec_tick   = 0:25:100;
set(gca, 'xtick', subframe_tick, 'xticklabel',msec_tick)
xlim([0 120])
title('Standard Basis')

xlabel('milliseconds');
orient landscape;
eval(sprintf('print -dpdf COB.pdf'))


