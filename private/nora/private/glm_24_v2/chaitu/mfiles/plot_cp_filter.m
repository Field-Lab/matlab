function out=plot_cp_filter(x,j,basepars,stimpars,trainpars)

Neff = size(trainpars.D,2);

out = 0;

cla;

if (Neff <= 1)
    
    if (basepars.nofilters_postspike == 0)
        return;
    end    
    xaxis = 0:trainpars.dt:(basepars.Mhist-1)*trainpars.dt;
else
    if (basepars.nofilters_coupling == 0)
        return;
    end
    xaxis = 0:trainpars.dt:(basepars.Mcoup-1)*trainpars.dt;
end


if (Neff <= 1)
    [~, PS] = get_model_filters(x,j,basepars,size(trainpars.D,2));
    % Plot the PS on short timescale
    plot(xaxis,PS), xlim([0 0.1]);
    title('Postspike filter: short timescale');
    return;
end


cpidx = get_pars_idx(basepars,j,Neff,'cp');
cp = basepars.coupling_basis*reshape(x(cpidx),basepars.nofilters_coupling,(size(trainpars.D,2)-1));

plot(xaxis,cp,'-','LineWidth',2);

scl = max(abs(cp(:)))/max(abs(basepars.coupling_basis(:)));

% Superimpose basis fns for illustrative purposes
hold on, plot(xaxis,scl*basepars.coupling_basis,'--');

%set(gca,'XScale','log');