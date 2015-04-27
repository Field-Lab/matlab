% print one-page summary of GLM params at the end of single_neuron.m (or glm_fit.m)
function printglmfit_AH2(Basepars,filename,savedir)
p0    = Basepars.p0;
p_opt = Basepars.p_opt;
dt    = Basepars.tstim / Basepars.spikebins_perstimframe;
%%
%%% Summary for non STA %%%
if strcmp(Basepars.k_filtermode, 'rk2')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    f_opt = opt_param.f;
    %    exitflag = opt_param.exitflag;
    %   figure(1), clf

        if ~Basepars.ps_FIX
        [K,PS] = get_model_filtersAH(p_opt,1,Basepars,1);
        else
            [K] = get_model_filtersAH(p_opt,1,Basepars,1);
            PS = p_opt(Basepars.paramind.PS)* (Basepars.ps_basis*Basepars.ps_Filter);
        end

        subplot(421)
        col_ax = [min(K(:)),max(K(:))];
        imagesc(K',col_ax)
        xlabel('Space [stixel]')
        ylabel('Time [frame]')
        title('Stimulus Filter (K)')
        if isfield(Basepars.paramind,'SR')
        	title(sprintf('Rect Scaled: %d',p_opt(Basepars.paramind.SR)));
        end
        
        
        
        
        colorbar

        subplot(422)
        [~,mx_idx] = max(abs(K(:)));
        [~,~,mx_fr] = ind2sub([11,11,20],mx_idx);
        imagesc(reshape(K(:,mx_fr),11,11),col_ax), colormap jet
        axis image
        set(gca,'xtick',5:5:11);
        set(gca,'ytick',5:5:11);
        title(sprintf('Stimulus filter at the peak (%dth frame)',mx_fr))
        colorbar

        LW = 2;
       
        subplot(423)
        cla
        t_ps = 0:dt:(Basepars.ps_timebins-1)*dt;
        t_ps = 1000*t_ps;
        plot(t_ps,exp(PS),'linewidth',LW)
        hold on
        plot([0,t_ps(end)],[1,1],'k--')
        ylabel('Gain')
        xlabel('Time [msec]')
        title('Post-Spike Filter (entire length)')
        xlim([0,t_ps(end)])

        subplot(424)
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
        xlim([0,25]); ylim([0,1.25]);
        title('Post-spike filter (0-15 msec)')

    %%%%%%% 
        if Basepars.Coupling  && length(Basepars.cp_Neighbors) > 0
            
            t_cp = 0:dt:(Basepars.cp_timebins-1)*dt;
            if Basepars.BiDirect_CP
                t_cp_negative = -(fliplr(t_cp)) -dt;
                t_cp          = [t_cp_negative, t_cp];
            end
            t_cp_millisecs = 1000* t_cp;
            
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
            
            
             zeroline_y_ON = linspace( min(exp(CP_ON(:))), max(exp(CP_ON(:))), 100);
             zeroline_y_OFF = linspace( min(exp(CP_ON(:))), max(exp(CP_OFF(:))), 100);
             zeroline_x = zeros(size(zeroline_y_ON));
             LW = 1.3;
             subplot(4,2,5); hold on; plot(t_cp_millisecs,exp(CP_ON),'linewidth',LW);   
                ylabel('Gain');
                xlabel('Time [msec]');
                title('On Parasol Filters');
                plot(zeroline_x, zeroline_y_ON, 'k');
                ylim([min(.5,min(exp(CP_ON(:)))), max(1.5,max(exp(CP_ON(:))))] ); 
               
                
            subplot(4,2,6); hold on; plot(t_cp_millisecs,exp(CP_OFF),'linewidth',LW);
                ylabel('Gain');
                xlabel('Time [msec]');
                title('Off Parasol Filters');
                plot(zeroline_x, zeroline_y_OFF, 'k');
                ylim([min(.5,min(exp(CP_OFF(:)))), max(1.5,max(exp(CP_OFF(:))))] ); 
            
                
             a = length(t_cp_millisecs);
             frontstart = round(.4*a);
             backend    = a -frontstart;
             fs = frontstart;
             be = backend;
                subplot(4,2,7); hold on; plot(t_cp_millisecs(fs:be),exp(CP_ON(fs:be,:)),'linewidth',LW);
                xlim ([t_cp_millisecs(fs) , t_cp_millisecs(be) ]);
                ylabel('Gain');
                xlabel('Time [msec]');
              %  title('ON Parasol CP');
                plot(zeroline_x, zeroline_y_ON, 'k');
                ylim([min(.5,min(exp(CP_ON(:)))), max(1.5,max(exp(CP_ON(:))))] ); 
               
                
            subplot(4,2,8); hold on; plot(t_cp_millisecs(fs:be),exp(CP_OFF(fs:be,:)),'linewidth',LW);
                xlim ([t_cp_millisecs(fs) , t_cp_millisecs(be) ]);
                ylabel('Gain');
                xlabel('Time [msec]');
              %  title('OFF Parasol CP');
                plot(zeroline_x, zeroline_y_OFF, 'k');
                ylim([min(.5,min(exp(CP_OFF(:)))), max(1.5,max(exp(CP_OFF(:))))] ); 
             
             
             
        end
%%%%%%%%%%%


        orient tall
        eval(sprintf('print -dpdf %s/%s_Rk2_%s_Summary.pdf',savedir, Basepars.fit_type,filename))
       
end


%%

if strcmp(Basepars.k_filtermode, 'STA');
    figure;
    
    
    K = p_opt(Basepars.paramind.L) * (Basepars.STA);
    if isfield(Basepars.paramind, 'FR')
        
        K = reshape(p_opt(Basepars.paramind.FR), 11,11);
    end
        
    
    
	subplot(421);
	col_ax = [min(K(:)),max(K(:))];
	imagesc(K',col_ax)
%	xlabel('Space [stixel]')
	ylabel('Time [frame]')
	title('STA Stimulus Filter (K)')

   if isfield(Basepars.paramind, 'FR')
        
       title('Nonlinearities')
    end
    colorbar
    
	subplot(422);
	c = 0;
	c = c+1;
	text(0, 1-0.1*c,sprintf('Baseline: %1.1e',p_opt(Basepars.paramind.MU)));
	c = c+3;
	text(0, 1-0.1*c,sprintf('STA Scaled: %1.1e',p_opt(Basepars.paramind.L)));
    c = c+5;
    
    if isfield(Basepars.paramind,'SR')
        	text(0, 1-0.1*c,sprintf('Rect Scaled: %1.1e',p_opt(Basepars.paramind.SR)));
    end
    axis off


  


    LW = 2;
    PS =  (Basepars.ps_basis)* p_opt(Basepars.paramind.PS);
    subplot(423)
	cla
	t_ps = 0:dt:(Basepars.ps_timebins-1)*dt;
	t_ps = 1000*t_ps;
	plot(t_ps,exp(PS),'linewidth',LW)
	hold on
	plot([0,t_ps(end)],[1,1],'k--')
	ylabel('Gain')
	xlabel('Time [msec]')
	title('Post-Spike Filter (entire length)')
	xlim([0,t_ps(end)])

     K = p_opt(Basepars.paramind.L) * (Basepars.STA);
	subplot(424)
	    %subplot(422)
        [~,mx_idx] = max(abs(K(:)));
        [~,~,mx_fr] = ind2sub([11,11,20],mx_idx);
        imagesc(reshape(K(:,mx_fr),11,11),col_ax), colormap jet
        axis image
        set(gca,'xtick',5:5:11);
        set(gca,'ytick',5:5:11);
        title(sprintf('Linear filter at the peak (%dth frame)',mx_fr))
        colorbar
    
	if Basepars.Coupling  && length(Basepars.cp_Neighbors) > 0
            
            t_cp = 0:dt:(Basepars.cp_timebins-1)*dt;
            if Basepars.BiDirect_CP
                t_cp_negative = -(fliplr(t_cp)) -dt;
                t_cp          = [t_cp_negative, t_cp];
            end
            t_cp_millisecs = 1000* t_cp;
            
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
            
            
             zeroline_y_ON = linspace( min(exp(CP_ON(:))), max(exp(CP_ON(:))), 100);
             zeroline_y_OFF = linspace( min(exp(CP_ON(:))), max(exp(CP_OFF(:))), 100);
             zeroline_x = zeros(size(zeroline_y_ON));
             LW = 1.3;
             subplot(4,2,5); hold on; plot(t_cp_millisecs,exp(CP_ON),'linewidth',LW);   
                ylabel('Gain');
                xlabel('Time [msec]');
                title('On Parasol Filters');
                plot(zeroline_x, zeroline_y_ON, 'k');
                ylim([min(.5,min(exp(CP_ON(:)))), max(1.5,max(exp(CP_ON(:))))] ); 
               
                
            subplot(4,2,6); hold on; plot(t_cp_millisecs,exp(CP_OFF),'linewidth',LW);
                ylabel('Gain');
                xlabel('Time [msec]');
                title('Off Parasol Filters');
                plot(zeroline_x, zeroline_y_OFF, 'k');
                ylim([min(.5,min(exp(CP_OFF(:)))), max(1.5,max(exp(CP_OFF(:))))] ); 
            
                
             a = length(t_cp_millisecs);
             frontstart = round(.4*a);
             backend    = a -frontstart;
             fs = frontstart;
             be = backend;
                subplot(4,2,7); hold on; plot(t_cp_millisecs(fs:be),exp(CP_ON(fs:be,:)),'linewidth',LW);
                xlim ([t_cp_millisecs(fs) , t_cp_millisecs(be) ]);
                ylabel('Gain');
                xlabel('Time [msec]');
              %  title('ON Parasol CP');
                plot(zeroline_x, zeroline_y_ON, 'k');
                ylim([min(.5,min(exp(CP_ON(:)))), max(1.5,max(exp(CP_ON(:))))] ); 
               
                
            subplot(4,2,8); hold on; plot(t_cp_millisecs(fs:be),exp(CP_OFF(fs:be,:)),'linewidth',LW);
                xlim ([t_cp_millisecs(fs) , t_cp_millisecs(be) ]);
                ylabel('Gain');
                xlabel('Time [msec]');
              %  title('OFF Parasol CP');
                plot(zeroline_x, zeroline_y_OFF, 'k');
                ylim([min(.5,min(exp(CP_OFF(:)))), max(1.5,max(exp(CP_OFF(:))))] ); 
             
             
             
        end
    orient tall
    eval(sprintf('print -dpdf %s/%s_STA_%s_Summary.pdf',savedir, Basepars.fit_type,filename))
    
   % eval(sprintf('print -dpdf %sScale_BiasEvolve.pdf',filename))
end
    
    
    
%%

    %%%%%%% 
detailedfilter_dir = sprintf('%s/FilterDetails',savedir);
if ~isdir(detailedfilter_dir), mkdir(detailedfilter_dir); end
    

%%%%%%%%%
%%%%%% Frame by Frame view %%%%%%%%%%
%{
figure;
if ~strcmp(Basepars.k_filtermode , 'STA')
    if ~Basepars.ps_FIX
        [K,PS] = get_model_filtersAH(p_opt,1,Basepars,1);
        %  [K0,PS0] = get_model_filtersAH(Basepars.p0,1,Basepars,1);
    else
        [K] = get_model_filtersAH(p_opt,1,Basepars,1);
        %[K0] = get_model_filtersAH(Basepars.p0,1,Basepars,1);
    end

    STA = Basepars.STA;
    
    [a,b] = max(abs(STA(:)));
    offon = sign(STA(b));  %%% determines if we have an on / off cell
    jitter_STA = std(STA(:));
    jitter_K = std(K(:));
    for plot_type = 1:3
        if plot_type == 1
            focus = 'FullFilter';
            STAplot = STA; Kplot = K;
        end
        
        if plot_type == 2  %%% plot onl positive parts
            STAplot = STA;  Kplot = K;
            negSTA = find(STA(:) <= jitter_STA);
            STAplot(negSTA) = 0;
            negK = find(K(:) <= jitter_K);
            Kplot(negSTA) = 0;
            if offon == 1
                focus = 'Center';
            else
                focus = 'Surround';
            end
        end 
        
        if plot_type == 3  %%% plot onl positive parts
            STAplot = STA;  Kplot = K;
            posSTA = find(STA(:) >= jitter_STA);
            STAplot(posSTA) = 0;
            posK = find(K(:) >= jitter_K);
            Kplot(posSTA) = 0;
            if offon == 1
                focus = 'Surround';
            else
                focus = 'Center';
            end
        end 
            
            
            
        for i_shift = 2:6
            shift = 3*(i_shift-1)+1;
            clf;
            col_ax = [ min(Kplot(:)) , max(Kplot(:)) ];
            subplot(321)
            imagesc(reshape(Kplot(:,shift),11,11),col_ax), colormap jet
            title(sprintf('K  %dth frame',shift)); colorbar

            subplot(323)
            imagesc(reshape(Kplot(:,shift+1),11,11),col_ax), colormap jet
            title(sprintf('K  %dth frame',shift+1)); colorbar

            subplot(325)
            imagesc(reshape(Kplot(:,shift+2),11,11),col_ax), colormap jet
            title(sprintf('K %dth frame',shift+2)); colorbar


            col_ax = [ min(STAplot(:)) , max(STAplot(:)) ];

            subplot(322)
            imagesc(reshape(STAplot(:,shift),11,11),col_ax), colormap jet
            title(sprintf('STA  %dth frame',shift)); colorbar

            subplot(324)
            imagesc(reshape(STAplot(:,shift+1),11,11),col_ax), colormap jet
            title(sprintf('STA  %dth frame',shift+1)); colorbar

            subplot(326)
            imagesc(reshape(STAplot(:,shift+2),11,11),col_ax), colormap jet
            title(sprintf('STA %dth frame',shift+2)); colorbar


            orient tall
            eval(sprintf('print -dpdf %s/%s_rk2_%s_%s%d.pdf',detailedfilter_dir, Basepars.fit_type,filename,focus,i_shift))

            %{
           % switch i
            %    case 1
                 eval(sprintf('print -dpdf %s_Filter%d.pdf',filename,type))
                     name = sprintf('%s_Filter%d',filename,type);
                     %saveas(gcf,name,'fig');
                     clear name;
                case 2
                 eval(sprintf('print -dpdf %s_FilterCenter%d.pdf',filename,type))
                    name = sprintf('%s_FilterCenter%d',filename,type);
                    %saveas(gcf,name,'fig');
                    clear name;
                case 3
                    name = sprintf('%s_FilterSurround%d',filename,type);
                    %saveas(gcf,name,'fig');
                    clear name;
                 eval(sprintf('print -dpdf %s_FilterSurround%d.pdf',filename,type))
            end
        end
            %}
        end
    end
end
%}


