function plot_basisfns(ktime_basis, kspace_basis, ps_basis, cp_basis, k_dt, ps_dt,xmode)

[mk nkb] = size(ktime_basis)
[mh nps] = size(ps_basis);
[mc nc] = size(cp_basis);


nplots = ceil(sqrt(1 + size(kspace_basis,2) + 3));


if (~exist('xmode','var'))
    xmode = 'log';
end

if (~isempty(ktime_basis))
    subplot(nplots,nplots,1);
    xaxis = 0:k_dt:(mk-1)*k_dt;
    plot(xaxis,ktime_basis), title('temporal linear filter basis fns'), xlabel('time (s)');
    set(gca,'XScale',xmode);
end

offset = 2;

if (~isempty(kspace_basis))
    for i=1:size(kspace_basis,2)
        subplot(nplots,nplots,offset);
        imagesc(reshape(kspace_basis(:,i),sqrt(size(kspace_basis,1)),sqrt(size(kspace_basis,1)))), axis image;
        title('spatial linear filter basis fns'), xlabel('time (s)');
        offset = offset + 1;
    end
end


if (~isempty(ps_basis))
    subplot(nplots,nplots,offset);
    xaxis= (0:ps_dt:(mh-1)*ps_dt);
    plot(xaxis,ps_basis), title('postspike filter basis fns'), xlabel('time (s)');
    set(gca,'XScale',xmode);
    offset = offset + 1;
end

if (~isempty(cp_basis))
    subplot(nplots,nplots,offset);
    xaxis = (0:ps_dt:(mc-1)*ps_dt);
    plot(xaxis,cp_basis), title('coupling filter basis fns'); xlabel('time (s)');
    set(gca,'XScale',xmode);
    offset = offset + 1;
end

if (~isempty(ps_basis) && ~isempty(ktime_basis))
    subplot(nplots,nplots,offset);
    plot(0:k_dt:(mk-1)*k_dt,ktime_basis), hold on, plot(0:ps_dt:(mh-1)*ps_dt,ps_basis,'--'), title('superimposed basis fns');
end