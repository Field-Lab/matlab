function out = plot_filter_prog(x,j,basepars,stimpars,trainpars)

out = 0;
%subplot(1,2,1);
%plot(x(end-9:end));


1;

kidx = get_pars_idx(basepars,j,size(trainpars.D,2),'k');

K = reshape(x(kidx),basepars.n,basepars.Mk);
cla;
if(basepars.n > 1)
    imagesc(K',[-1 1].*max(eps,max(abs(K(:))))), box on, colormap gray, axis tight
    xlabel('Space [stixel]')
    ylabel('Time [frame]')
else
    plot(0:stimpars.dt:(basepars.Mk-1)*stimpars.dt,K,'.-');
end
title(sprintf('b=%0.5f stimulus filter norm: %0.3f/%0.3f',x(get_pars_idx(basepars,j,size(trainpars.D,2),'b')),norm(K,'fro')*sqrt(stimpars.dt),sum(abs(K(:)))*stimpars.dt));

%subplot(1,2,2);
%plot(x(end-9:end));
%pause(0.1);

