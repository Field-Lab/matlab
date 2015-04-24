function out=plot_ps_filter(x,j,basepars,stimpars,trainpars)

out = 0;
cla;
xaxis = 0:trainpars.dt:(basepars.Mhist-1)*trainpars.dt;


offset = get_pars_idx(basepars,j,size(trainpars.D,2),'ps');
if (isempty(offset))
    return;
end
offset = offset(1)-1;

psidx = offset+1:offset+basepars.nofilters_postspike;
ps = basepars.postspike_basis*x(psidx);

plot(xaxis,ps,'LineWidth',2);

scl = max(abs(ps))/max(abs(basepars.postspike_basis(:)));

% Superimpose basis fns for illustrative purposes
hold on, plot(xaxis,scl*basepars.postspike_basis,'--');

1;
if 0
cutoff = find_last_nonzero(ps,10^(-3));
%cutoff
if (~isempty(cutoff) && cutoff > 1 && cutoff <= length(xaxis))
    xlim([0 xaxis(cutoff)]);
end
end
ylim([-0.1 0.1]);
%set(gca,'XScale','log');
title(sprintf('Postspike filter: long timescale. Area=%0.3f',sum(ps)*trainpars.dt));