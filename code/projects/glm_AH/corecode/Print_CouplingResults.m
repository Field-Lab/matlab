%%% More Model Analysis
%%% first load up the rk2_1_cpOFF  blahblah
fn_save = Basepars.fn_save;
 
p0     = opt_param.p0;
p_opt  = opt_param.p; 

 
n_off = length ( Basepars.OFF_Neighbor);
n_on  = length ( Basepars.ON_Neighbor);
 

CP_OFF = zeros( max(size(Basepars.cp_basis)) , n_off );
CP_ON  = zeros( max(size(Basepars.cp_basis)) , n_on);
 
CP_OFF0 = zeros( max(size(Basepars.cp_basis)) , n_off );
CP_ON0  = zeros( max(size(Basepars.cp_basis)) , n_on);
 
for i = 1 : n_off
    cid        = Basepars.OFF_Neighbor(i) ;
    n_idx      = find(Basepars.cp_Neighbors == cid);
    paramshift = min(Basepars.paramind.CP) ;
    
    %%%%%%%%%%%% 
    paramind     = [ (paramshift + (n_idx-1)*Basepars.cp_filternumber) : (paramshift + n_idx*Basepars.cp_filternumber -1) ];
    CP_OFF(:,i)  = Basepars.cp_basis * (p_opt(paramind ));
    CP_OFF0(:,i) = Basepars.cp_basis * (p0(paramind));
end
for i = 1 : n_on
    cid        = Basepars.ON_Neighbor(i) ;
    n_idx      = find(Basepars.cp_Neighbors == cid);
    paramshift = min(Basepars.paramind.CP) ;
    
    %%%%%%%%%%%% 
    paramind     = [ (paramshift + (n_idx-1)*Basepars.cp_filternumber) : (paramshift + n_idx*Basepars.cp_filternumber -1) ];
    CP_ON(:,i)  = Basepars.cp_basis * (p_opt(paramind ));
    CP_ON0(:,i) = Basepars.cp_basis * (p0(paramind));
end
     
 
figure;  LW = 2;
subplot(2,2,1); hold on; plot(exp(CP_ON),'linewidth',LW);   
    ylabel('Gain');
    xlabel('Time [msec]');
    title('ON Parasol Opt');
subplot(2,2,2); hold on; plot(exp(CP_OFF),'linewidth',LW);
    ylabel('Gain');
    xlabel('Time [msec]');
    title('OFF_Parasol Opt');
subplot(2,2,3); hold on; plot(exp(CP_ON0),'linewidth',LW);   
    ylabel('Gain');
    xlabel('Time [msec]');
    title('ON Parasol Init');
    
subplot(2,2,4); hold on; plot(exp(CP_OFF0),'linewidth',LW);
    ylabel('Gain');
    xlabel('Time [msec]');
    title('OFF Parasol Init');
eval(sprintf('print -dpdf %s_Coupling.pdf',fn_save))
  
  %{
for i = 2:3
    switch i
        case 1    
            K  = Model.params.LinearFilter;
            PS = Basepars.ps_basis * Model.params.ps;
            if Basepars.Coupling
                shift = min (cp_idx{1}) - 1; 
                for iNeigh = 1 : nNeighbors
                    CP(:,iNeigh) = Basepars.cp_basis *  (Model.params.cp ( cp_idx{iNeigh} - shift ) );
                end
            end
        case 2
            K  = p0(s1_idx)*p0(t1_idx)' + p0(s2_idx)*p0(t2_idx)';
            PS = Basepars.ps_basis * p0(ps_idx)
            if Basepars.Coupling
                shift = 0 ;
                for iNeigh = 1 : nNeighbors
                    CP(:,iNeigh) = Basepars.cp_basis *  ( p0( cp_idx{iNeigh} - shift ) );
                end
            end
        case 3
            K  = p_opt(s1_idx)*p_opt(t1_idx)' + p_opt(s2_idx)*p_opt(t2_idx)';
            PS =  Basepars.ps_basis*p_opt(ps_idx);
            if Basepars.Coupling
                shift = 0 ;
                for iNeigh = 1 : nNeighbors
                    CP(:,iNeigh) = Basepars.cp_basis *  ( p_opt( cp_idx{iNeigh} - shift ) );
                end
            end
    end
            
    figure;
    subplot(221)
    col_ax = [min(K(:)),max(K(:))];
    imagesc(K',col_ax)
    xlabel('Space [stixel]')
    ylabel('Time [frame]')
    title('Stimulus filter (K)')
    colorbar

    subplot(223)
    [~,mx_idx] = max(abs(K(:)));
    [~,~,mx_fr] = ind2sub([15,15,30],mx_idx);
    imagesc(reshape(K(:,mx_fr),15,15),col_ax), colormap jet
    axis image
    set(gca,'xtick',5:5:15);
    set(gca,'ytick',5:5:15);
    title(sprintf('Stimulus filter at the peak (%dth frame)',mx_fr))
    colorbar

    LW = 2;

    subplot(222)
    cla
    t_ps = 0:Trainpars.dt:(Basepars.ps_timebins-1)*Trainpars.dt;
    t_ps = 1000*t_ps;
    plot(t_ps,exp(PS),'linewidth',LW)
    hold on
    plot([0,t_ps(end)],[1,1],'k--')
    ylabel('Gain')
    xlabel('Time [msec]')
    title('Post-spike filter (entire length)')
    xlim([0,t_ps(end)])
    if Basepars.Coupling
        subplot(224)
        plot(CP, 'linewidth',LW)
    end
  %  plot([0,t_ps(end)],[1,1],'k--')
  %  hold on
  %  for k = 1:Basepars.ps_filternumber
  %     psc = Basepars.ps_basis(:,k)*p_opt(1+Basepars.k_spacepixels+k);
  %    plot(t_ps,exp(psc),'k','linewidth',LW/3)
  %  end
  %  plot(t_ps,exp(PS),'linewidth',LW)
  %  ylabel('Gain')
  %  xlabel('Time [msec]')
  %  xlim([0,15])
  %  title('Post-spike filter (0-15 msec)')  
    switch i
        case 1
            eval(sprintf('print -dpdf %s_SimulationFilter.pdf',fn_save))
        case 2
            eval(sprintf('print -dpdf %s_InitialGuess.pdf',fn_save))
        case 3
            eval(sprintf('print -dpdf %s_OptimizedGuess.pdf',fn_save))
    end
end
%}