function out = plot_sep_filter_prog(x,j,basepars,stimpars,trainpars)

out = 0;
%subplot(1,2,1);
%plot(x(end-9:end));

neff = size(trainpars.D,2);
offset = get_pars_idx(basepars,j,neff,'k');
offset = offset(1)-1;
[space1_idx temp1_idx space2_idx temp2_idx] = get_sep_filt_idces(offset,basepars);

switch(basepars.filtermode)
    case {'sep_raw','rk2'}
        K = x(space1_idx)*x(temp1_idx)' + x(space2_idx)*x(temp2_idx)';
    case 'sep_basis'
        K = (basepars.kspace_basis*x(space1_idx))*(basepars.ktime_basis*x(temp1_idx))' + (basepars.kspace_basis*x(space2_idx))*(basepars.ktime_basis*x(temp2_idx))';
end
cla;
1;

%imgb(K,0,0,[sqrt(basepars.n) sqrt(basepars.n)])
imagesc(K', [-1 1].*max(eps,max(abs(K(:))))), axis tight off; colormap gray, xlabel('Space [stixel]'),ylabel('Time [frame]')
title(sprintf('b=%0.5f,stim norm=%0.5f',x(get_pars_idx(basepars,j,neff,'b')), norm(K,'fro')*sqrt(stimpars.dt)));
%showIm(K);%, colorbar;
%pause(1.0);
%subplot(1,2,2);
%plot(x(end-9:end));
%pause;

