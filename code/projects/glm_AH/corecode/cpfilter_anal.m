%%% Just some crappy diagnostics 
%%% JUST LOAD THE Model_1_ps12.mat   


p0 = opt_param.p0;
p_opt = opt_param.p;


msecdur = 80;  %%% millisecs.. must be less than 100
bindur  = round( Basepars.spikebins_perstimframe* 120*(msecdur /1000) );



%figure; subplot(2,1,1); plot(Basepars.ps_basis*PSpars);
%subplot(2,1,2); plot(Basepars.ps_basis*orig_PSPars);


Pvec_progression = [Track_Progress.p];
iterations = Track_Progress(1).counter;

if ~Basepars.ps_FIX
    PSpars = opt_param.p(Basepars.paramind.PS);
    orig_PSPars = Basepars.p0(Basepars.paramind.PS);
    multiple = floor(iterations/6);
    figure;
    for i = 1 : 6
        subplot(3,2,i);
        ind = 1 + 6 *i;
        if ind > iterations
            ind = iterations;
        end
        psfilter = (Basepars.ps_basis*Pvec_progression(Basepars.paramind.PS, ind));
        l = length(psfilter);
        plot(psfilter(1:floor(l/4)));
    end
    psname = sprintf('%s_PSconv',Basepars.fn_save);
    eval(sprintf('print -dpdf %s.pdf',psname))
end


    figure; subplot(2,1,1);
    plot(Pvec_progression(1,:)); xlabel('mu');
    subplot(2,1,2);
    plot(Pvec_progression(2,:)); xlabel('scaleSTA');
    mu_sta_name = sprintf('%s_Mu_STA',Basepars.fn_save);
    eval(sprintf('print -dpdf %s.pdf',mu_sta_name))
  
%%%%
if Basepars.Coupling && length(Basepars.cp_Neighbors) > 0 
    clear cp_perNeighbor
    nNeighbors = length(Basepars.cp_Neighbors);
    cpname = sprintf('cell%d_binperframe%d_cpbasis%d_fitsec%d' ,Basepars.headid,Basepars.spikebins_perstimframe,...
        Basepars.cp_filternumber,round(Basepars.maxt / 120));

    cp_ind = Basepars.paramind.CP;
    cp_optparams = p_opt(cp_ind);


    %%%%%  THIS IS FOR ORIGINAL TO OTHER COMPARISON
    %{   

    cpfilters = zeros(nNeighbors , max(size(Basepars.cp_basis)));
    for i = 1:nNeighbors
        cp_perNeighbor{i}.params = cp_optparams((i-1)*Basepars.cp_filternumber + 1 : i *Basepars.cp_filternumber);
        cp_perNeighbor{i}.cpfilter = Basepars.cp_basis * cp_perNeighbor{i}.params;
        cpfilters(i,:) = cp_perNeighbor{i}.cpfilter';
    end
    figure;
    subplot(2,1,1); hold on
    plot(exp(cpfilters(:,1:bindur))');xlabel('OptimumCPFiltersGain'); hold off

    clear i

    cp_origparams  = Basepars.p0(cp_ind);
    cp_origfilters = zeros(nNeighbors , max(size(Basepars.cp_basis)));
    for i = 1:nNeighbors
        cp_perNeighbor_params= cp_origparams((i-1)*Basepars.cp_filternumber + 1 : i *Basepars.cp_filternumber);
        cp_perNeighbor_cpfilter = Basepars.cp_basis * cp_perNeighbor_params;
        cp_origfilters(i,:) = cp_perNeighbor_cpfilter';
    end

    subplot(2,1,2);  hold on
    plot(exp(cp_origfilters(:,1:bindur))');  xlabel('OriginalCPFiltersGain'); hold off
    %}



    %%% OFF ON AND JUNK
    % OFF FIRST

    figure;  ymin = .7; ymax =1.8;
    offcp_filters = zeros(length(Basepars.OFF_Neighbor) , max(size(Basepars.cp_basis)));
    for i = 1:length(Basepars.OFF_Neighbor)
        nind = find(Basepars.cp_Neighbors == Basepars.OFF_Neighbor(i));

        cp_params          =  cp_optparams((nind-1)*Basepars.cp_filternumber + 1 : nind *Basepars.cp_filternumber);
        offcp_filters(i,:) = (Basepars.cp_basis * cp_params)';
    end
    offcp_filters = offcp_filters(:,1:bindur)'

    subplot(3,1,1); hold on; ylim([ymin,ymax]); 
    plot(exp(offcp_filters)); plot(mean(exp(offcp_filters),2),'k.');
    ylabel('OFF'); xlabel('gain'); hold off;

    oncp_filters = zeros(length(Basepars.ON_Neighbor) , max(size(Basepars.cp_basis)));
    for i = 1:length(Basepars.ON_Neighbor)
        nind = find(Basepars.cp_Neighbors == Basepars.ON_Neighbor(i));

        cp_params          =  cp_optparams((nind-1)*Basepars.cp_filternumber + 1 : nind *Basepars.cp_filternumber);
        oncp_filters(i,:) = (Basepars.cp_basis * cp_params)';
    end
    oncp_filters = oncp_filters(:,1:bindur)'
    
    subplot(3,1,2); hold on; ylim([ymin,ymax]); 
    plot(exp(oncp_filters));plot(mean(exp(oncp_filters),2),'k.');
    ylabel('ON'); xlabel('gain'); hold off;


    other_id = setdiff(Basepars.cp_Neighbors , union( Basepars.ON_Neighbor , Basepars.OFF_Neighbor) );
    other_filters = zeros(  length(other_id)  , max(size(Basepars.cp_basis))  );
    for i = 1 : length(other_id)
          nind = find(Basepars.cp_Neighbors == other_id(i) );
          cp_params          =  cp_optparams((nind-1)*Basepars.cp_filternumber + 1 : nind *Basepars.cp_filternumber);
          other_filters(i,:) = (Basepars.cp_basis * cp_params)';
    end
    other_filters = other_filters(:,1:bindur)'

    subplot(3,1,3); hold on; ylim([ymin,ymax]); 
    plot(exp(other_filters));
    ylabel('Distal'); xlabel('gain'); hold off;
    
    cp_name = sprintf('%s_CP',Basepars.fn_save);
    eval(sprintf('print -dpdf %s.pdf',cp_name))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
eval(sprintf('print -dpdf %s.pdf',cpname))
%%% other filters  %%%

% Linear Filter
k_ind = get_pars_idx2AH(Basepars,'kfilter');
figure; subplot(2,1,1)
subplot(2,1,1);hold on; plot(p_opt(k_ind));
xlabel('FinalKFilter'); hold off;
subplot(2,1,2);plot(p0(k_ind)); hold on;
xlabel('OriginalKfilter'); hold off;

kname = sprintf('CollapsetoRK1_cell%d_bin%d_cpon_fitsec%d' ,Basepars.headid,Basepars.spikebins_perstimframe,...
    round(Basepars.maxt / 120));
eval(sprintf('print -dpdf %s.pdf',kname))
%}