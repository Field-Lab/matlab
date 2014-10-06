% print one-page summary of GLM params at the end of single_neuron.m (or glm_fit.m)
function print_glm_fitAH(cid,opt_param,Basepars,dt,fn_save,savedir,Track_Progress)

if ~strcmp(Basepars.k_filtermode, 'STA')
    for type = 1 : 2
        if type == 1
            p_opt = opt_param.p;
        elseif type == 2
            p_opt = opt_param.p0;
        end
        f_opt = opt_param.f;
        exitflag = opt_param.exitflag;

        figure(1), clf

        if exitflag > 0
           opt_end = 'optimization was converged.';
        elseif exitflag == 0
           opt_end = sprintf('max iter was reached before the convergence.');
        elseif exitflag == -1
           opt_end = sprintf('optimization was terminated with some error.');
        end

       subplot(321)
       c = 0;
       c = c+1;
       text(0, 1-0.1*c,sprintf('Cell ID: %4d',cid))
       c = c+1;
       text(0, 1-0.1*c,sprintf('Baseline: %1.1e',p_opt(1)))
       c = c+1;
     %  text(0, 1-0.1*c,sprintf('Negative LL: %1.1e',f_opt))
     %  c = c+1;
       text(0, 1-0.1*c,sprintf('%s',opt_end))
       axis off

        %-- extract filters
        if ~Basepars.ps_FIX
        [K,PS] = get_model_filtersAH(p_opt,1,Basepars,1);
        else
            [K] = get_model_filtersAH(p_opt,1,Basepars,1);
            PS = p_opt(Basepars.paramind.PS)* (Basepars.ps_basis*Basepars.ps_Filter);
        end

        subplot(323)
        col_ax = [min(K(:)),max(K(:))];
        imagesc(K',col_ax)
        xlabel('Space [stixel]')
        ylabel('Time [frame]')
        title('Stimulus filter (K)')
        colorbar

        subplot(325)
        [~,mx_idx] = max(abs(K(:)));
        [~,~,mx_fr] = ind2sub([15,15,30],mx_idx);
        imagesc(reshape(K(:,mx_fr),15,15),col_ax), colormap jet
        axis image
        set(gca,'xtick',5:5:15);
        set(gca,'ytick',5:5:15);
        title(sprintf('Stimulus filter at the peak (%dth frame)',mx_fr))
        colorbar

        LW = 2;

        subplot(324)
        cla
        t_ps = 0:dt:(Basepars.ps_timebins-1)*dt;
        t_ps = 1000*t_ps;
        plot(t_ps,exp(PS),'linewidth',LW)
        hold on
        plot([0,t_ps(end)],[1,1],'k--')
        ylabel('Gain')
        xlabel('Time [msec]')
        title('Post-spike filter (entire length)')
        xlim([0,t_ps(end)])

        subplot(326)
        cla
        plot([0,t_ps(end)],[1,1],'k--')
        hold on
        for k = 1:Basepars.ps_filternumber
           psc = Basepars.ps_basis(:,k)*p_opt(1+Basepars.k_spacepixels+k);
           plot(t_ps,exp(psc),'k','linewidth',LW/3)
        end
        plot(t_ps,exp(PS),'linewidth',LW)
        ylabel('Gain')
        xlabel('Time [msec]')
        xlim([0,15])
        title('Post-spike filter (0-15 msec)')


        orient tall
        if type == 1
            %saveas(gcf , fn_save, 'fig');
            eval(sprintf('print -dpdf %s.pdf',fn_save))
        elseif type == 2
            name = sprintf('%sinitialguess',fn_save);
            %saveas(gcf,name,'fig');
            clear name;
            eval(sprintf('print -dpdf %s/%sinitialguess.pdf',savedir,fn_save))
        end

    end
end



if strcmp(Basepars.k_filtermode, 'STA');
    figure;
    
          f_opt = opt_param.f;
        exitflag = opt_param.exitflag;

        figure(1), clf

        if exitflag > 0
           opt_end = 'optimization was converged.';
        elseif exitflag == 0
           opt_end = sprintf('max iter was reached before the convergence.');
        elseif exitflag == -1
           opt_end = sprintf('optimization was terminated with some error.');
        end

       subplot(311)
       bias      = p_evolve (Basepars.paramind.MU,:);
    STAscale  = p_evolve(Basepars.paramind.L,:);
       c = 0;
       c = c+1;
       text(0, 1-0.1*c,sprintf('Cell ID: %d',cid))
       c = c+1;
       text(0, 1-0.1*c,sprintf('Baseline: %d',bias))
       c = c+1;
       text(0, 1-0.1*c,sprintf('STA_Scale: %d', STAscale))
       c = c+ 1;
       
       text(0, 1-0.1*c,sprintf('Negative LL: %1.1e',f_opt))
       
       c = c+1;
       text(0, 1-0.1*c,sprintf('%s',opt_end))
       axis off
    
    
    
    
%    subplot(3,1,1);
%    count = Track_Progress(1).counter;
    
 %   p_evolve = [Track_Progress(1:count).p];
    
    
    
  %  plot(bias,'b');
  %  hold on;
  %  plot(STAscale,'r');
  %  xlabel('iterations');
    subplot(3,1,2);
    imagesc(STAscale(end)* Basepars.STA'); colorbar
    
    if ~Basepars.ps_FIX
        subplot(3,1,3);
        p_opt   = opt_param.p;
        psparams =  p_opt(Basepars.paramind.PS);
        psfilter = Basepars.ps_basis * psparams;
        psgain = exp(psfilter);
        plot(psgain, 'm*');
    end
    
    
    if Basepars.ps_FIX
        subplot(3,2,5);
        plot(Basepars.ps_FIXconvolved);
        xlabel('shouldbezero')
        subplot(3,2,6);
        PS_scale = p_evolve(Basepars.paramind.PS,:);
        plot(PS_scale);
        xlabel('PS Scale')
    end
        
    

    
    eval(sprintf('print -dpdf %s/%sScale_BiasEvolve.pdf',savedir,fn_save))
end
    
    
    


    %%%%%%% 
    if Basepars.Coupling  && length(Basepars.cp_Neighbors) > 0
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


        name = sprintf('%s/%s_Coupling',savedir,fn_save);
        %saveas(gcf,name,'fig');
        clear name;
    end



%%%%%%%%%
%%

if ~strcmp(Basepars.k_filtermode , 'STA')
    if ~Basepars.ps_FIX
        [K,PS] = get_model_filtersAH(p_opt,1,Basepars,1);
      %  [K0,PS0] = get_model_filtersAH(Basepars.p0,1,Basepars,1);
        else
            [K] = get_model_filtersAH(p_opt,1,Basepars,1);
            %[K0] = get_model_filtersAH(Basepars.p0,1,Basepars,1);
    end
    STA = Basepars.STA;


    for i = 1:3

        [~,mx_idx] = max(abs(K(:)));
        [~,~,mx_fr] = ind2sub([15,15,30],mx_idx);

        if i == 3
            [~,mx_idx] = max((K(:)));
            [~,~,mx_fr] = ind2sub([15,15,30],mx_idx);
        end

        for type = 1:2
            figure;
            if type ==  1
                shift = -3;
            elseif type ==2
                shift = 0;
            end


            switch i
                case 1
                    col_ax = [ min(K(:)) , max(K(:)) ];
                case 2
                    col_ax = [ min(K(:)) ,  -std(K(:)) ];
                case 3
                    col_ax = [ std(K(:))         , max(K(:)) ];
            end

            subplot(321)
            imagesc(reshape(K(:,mx_fr-shift),15,15),col_ax), colormap jet
            title(sprintf('K  %dth frame',mx_fr-shift)); colorbar

            subplot(323)
            imagesc(reshape(K(:,mx_fr-shift+1),15,15),col_ax), colormap jet
            title(sprintf('K  %dth frame',mx_fr-shift+1)); colorbar

            subplot(325)
            imagesc(reshape(K(:,mx_fr-shift+2),15,15),col_ax), colormap jet
            title(sprintf('K %dth frame',mx_fr-shift+2)); colorbar


           switch i
                case 1
                    col_ax = [ min(STA(:)) , max(STA(:)) ];
                case 2
                    col_ax = [ min(STA(:)) ,        - std(STA(:)) ];
                case 3
                    col_ax = [ std(STA(:))           , max(STA(:)) ];
            end

            subplot(322)
            imagesc(reshape(STA(:,mx_fr-shift),15,15),col_ax), colormap jet
            title(sprintf('STA  %dth frame',mx_fr-shift)); colorbar

            subplot(324)
            imagesc(reshape(STA(:,mx_fr-shift+1),15,15),col_ax), colormap jet
            title(sprintf('STA  %dth frame',mx_fr-shift+1)); colorbar

            subplot(326)
            imagesc(reshape(STA(:,mx_fr-shift+2),15,15),col_ax), colormap jet
            title(sprintf('STA %dth frame',mx_fr-shift+2)); colorbar


            orient tall
            switch i
                case 1
                 eval(sprintf('print -dpdf %s_Filter%d.pdf',fn_save,type))
                     name = sprintf('%s_Filter%d',fn_save,type);
                     %saveas(gcf,name,'fig');
                     clear name;
                case 2
                 eval(sprintf('print -dpdf %s_FilterCenter%d.pdf',fn_save,type))
                    name = sprintf('%s_FilterCenter%d',fn_save,type);
                    %saveas(gcf,name,'fig');
                    clear name;
                case 3
                    name = sprintf('%s_FilterSurround%d',fn_save,type);
                    %saveas(gcf,name,'fig');
                    clear name;
                 eval(sprintf('print -dpdf %s_FilterSurround%d.pdf',fn_save,type))
            end
        end
    end
end


