% Plot the gradient values

function out = plot_gradient(x,optimValues,basepars,stimpars,trainpars)

1;

gradient = optimValues.gradient;


gradient = gradient(:);

1;
%if (isnan(optimValues.ratio))
%    out = 1;
%    return;
%end

out = 0;

Neff = size(trainpars.D,2);

if (isempty(gradient))
    return;
end

b_idx = get_pars_idx(basepars,1,Neff,'b');

switch (basepars.filtermode)

    case 'sep_basis'

        [nspace ntime] = get_nspacetime(basepars);
        offset = get_pars_idx(basepars,1,Neff,'k');
        offset = offset(1)-1;
        [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(offset,basepars);
        
    
    case 'sep_raw'

        [nspace ntime] = get_nspacetime(basepars);
        offset = get_pars_idx(basepars,1,Neff,'k');
        offset = offset(1)-1;
        [s1_idx t1_idx s2_idx t2_idx] = get_sep_filt_idces(offset,basepars);
        
    case 'nonsep'
        
        k_idx = get_pars_idx(basepars,1,Neff,'k');
        
    case 'fixfilt'
        kn_idx = get_pars_idx(basepars,1,Neff,'k');
        
end

ps_idx = get_pars_idx(basepars,1,Neff,'ps');

cp_idx = reshape(get_pars_idx(basepars,1,Neff,'cp'),Neff-1,basepars.nofilters_coupling);

1;
cla;
bar(1,log(abs(gradient(b_idx))),'b');

offset = 1;

switch(basepars.filtermode)
    
    case 'sep_basis'
        
        for j=1:nspace
            hold on, bar(offset+1,log(abs(gradient(s1_idx(j)))),'r');
            offset = offset + 1;
        end
        for j=1:ntime
            hold on, bar(offset+1,log(abs(gradient(t1_idx(j)))),'y');
            offset = offset + 1;
        end        
        for j=1:nspace
            hold on, bar(offset+1,log(abs(gradient(s2_idx(j)))),'r');
            offset = offset + 1;
        end        
        for j=1:ntime
            hold on, bar(offset+1,log(abs(gradient(t2_idx(j)))),'y');
            offset = offset + 1;
        end    
        
    case 'sep_raw'
        hold on, bar(offset+1,log(mean(abs(gradient(s1_idx)))),'r');
        offset = offset + 1;%basepars.n;
        hold on, bar(offset+1,log(mean(abs(gradient(t1_idx)))),'y');
        offset = offset + 1;%basepars.Mk;
        hold on, bar(offset+1,log(mean(abs(gradient(s2_idx)))),'r');
        offset = offset + 1;%basepars.n;
        hold on, bar(offset+1,log(mean(abs(gradient(t2_idx)))),'y');
        offset = offset + 1;%basepars.Mk;        
    case 'nonsep'
        
        hold on, bar(offset+1:offset+1,log(mean(abs(gradient(k_idx)))),'r');
        offset = offset + 1;
        
    case 'fixfilt'
        hold on, bar(offset+1:offset+1,log(abs(gradient(kn_idx))),'r');
        offset = offset + 1;
        
end
if (basepars.nofilters_postspike > 0)
    hold on, bar(offset+1:offset+length(ps_idx),log(abs(gradient(ps_idx))),'g');
    offset = offset+length(ps_idx);
end

if (Neff>1 && basepars.nofitlers_coupling > 0)
    for j=1:basepars.nofilters_coupling
        hold on, bar(offset+1,log(mean(abs(gradient(cp_idx(j,:))))),'c');
        offset = offset + 1;
    end
end