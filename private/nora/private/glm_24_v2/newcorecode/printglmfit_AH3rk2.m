% TRYING TO MAKE MORE MODULAR
% Printglmfit_AH2 is good and fast and works .. just not modular

function printglmfit_AH3rk2(Basepars,savedir,filename)
%{
opt_param.p     = p_opt;
opt_param.p0    = p0;
opt_param.f     = f_opt;
opt_param.g     = g_opt;
opt_param.H     = H_opt;
opt_param.exitflag = eflag;
opt_param.output   = output;
opt_param.lcifs = lcifs;
% trying to get here .. not there yet;
printglmfit_AH3rk2(cid, opt_param,paramind,ps_basis, cp_basis, tstim, dt, filename)
%}

%{
bins = Basepars.maxt * 10;
avglogprob = opt_param.f / bins;
%}
opt_param = Basepars.opt_param;
info.cid = Basepars.headid;
info.exp_nm = Basepars.exp_nm;
info.fittype = Basepars.fit_type;
info.ctype = Basepars.celltype;
info.conemodel = Basepars.conemodel;


    

Z = Basepars.paramind;
p_opt = opt_param.p;
ps_basis = Basepars.ps_basis;
dt = Basepars.dt;
tstim =Basepars.tstim;
cid = Basepars.headid;
opt_param = Basepars.opt_param;
%clear Trainpars cid eflag output
%%
%%% Summary for non STA %%%
figure;


if strcmp(Basepars.k_filtermode, 'fixedSP')
    MU = p_opt(Z.MU);
    K  = (Basepars.spfilter)*(p_opt(Z.L)');
    PS = ps_basis * p_opt(Z.PS);
    clf
    subplot(3,2,1)
    set(gca, 'fontsize', 12)
    c = 0;
    c = c+1;
    text(0, 1-0.1*c,sprintf('%s Fit by %s',info.exp_nm, info.fittype))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Cone Model: %s',info.conemodel))
    c = c+1;
    text(0, 1-0.1*c,sprintf('%s  Cell ID: %4d', info.ctype,info.cid))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Tonic drive with gray stim: %d hz',round(exp(MU)))  )
    c = c+1;
	text(0, 1-0.1*c,sprintf('Optimum fmin: %1.5e',opt_param.f))
	c = c+1;
    text(0, 1-0.1*c,sprintf('fminunc exit flag: %d',opt_param.exitflag))
    if isfield(Basepars.GLMPars, 'rect') && Basepars.GLMPars.rect
        c = c+1;
        text(0, 1-0.1*c,sprintf('Stim is Rectified with %s',Basepars.GLMPars.rect_type));
    end
	axis off

	subplot(3,2,2)  % label in 50 msec intervals
    spfilter = (Basepars.spfilter);
    klen = sqrt(length(spfilter));
    imagesc(reshape(spfilter,[klen,klen]));
    colorbar
    
    
    LW = 2;
	subplot(3,2,4)   
    set(gca, 'fontsize', 12)
    timebins = length(Z.L);
    t_tf = 1000*(0:dt:(timebins-1)*dt);
    plot(t_tf,p_opt(Z.L),'linewidth',LW)
	hold on
	plot([0,t_tf(end)],[0,0],'k--') % black reference line 
    xlabel('Time [msec]')
    xlim([0,t_tf(end)])
	title(sprintf('Time cours .. temporal filter'));
    
    

       
    subplot(323)
    set(gca, 'fontsize', 12)
    timebins = length(PS);
    t_ps = 1000*(0:dt:(timebins-1)*dt);
	plot(t_ps,exp(PS),'linewidth',LW)
	hold on
	plot([0,t_ps(end)],[1,1],'k--') % black reference line 
    ylabel('Gain')
    xlabel('Time [msec]')
    title('Post-Spike Filter')
    xlim([0,t_ps(end)])

        

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
           % CP_OFF0 = zeros( max(size(Basepars.cp_basis)) , n_off );
           % CP_ON0  = zeros( max(size(Basepars.cp_basis)) , n_on);

            for i = 1 : n_off
                cid        = Basepars.OFF_Neighbor(i) ;
                n_idx      = find(Basepars.cp_Neighbors == cid);
                paramshift = min(Basepars.paramind.CP) ;

                %%%%%%%%%%%% 
                paramind     = [ (paramshift + (n_idx-1)*Basepars.cp_filternumber) : (paramshift + n_idx*Basepars.cp_filternumber -1) ];
                CP_OFF(:,i)  = Basepars.cp_basis * (p_opt(paramind ));
             %   CP_OFF0(:,i) = Basepars.cp_basis * (p0(paramind));
            end
            for i = 1 : n_on
                cid        = Basepars.ON_Neighbor(i) ;
                n_idx      = find(Basepars.cp_Neighbors == cid);
                paramshift = min(Basepars.paramind.CP) ;

                %%%%%%%%%%%% 
                paramind     = [ (paramshift + (n_idx-1)*Basepars.cp_filternumber) : (paramshift + n_idx*Basepars.cp_filternumber -1) ];
                CP_ON(:,i)  = Basepars.cp_basis * (p_opt(paramind ));
            %    CP_ON0(:,i) = Basepars.cp_basis * (p0(paramind));
            end
            
            
             zeroline_y_ON = linspace( min(exp(CP_ON(:))), max(exp(CP_ON(:))), 100);
             zeroline_y_OFF = linspace( min(exp(CP_ON(:))), max(exp(CP_OFF(:))), 100);
             zeroline_x = zeros(size(zeroline_y_ON));
             LW = 1.3;
             subplot(3,2,5); set(gca, 'fontsize', 12);
             hold on; plot(t_cp_millisecs,exp(CP_ON),'linewidth',LW);   
                ylabel('Gain');
                xlabel('Time [msec]');
                title('On Parasol Filters');
                plot(zeroline_x, zeroline_y_ON, 'k');
                ylim([min(.5,min(exp(CP_ON(:)))), max(1.5,max(exp(CP_ON(:))))] ); 
               
                
            subplot(3,2,6); set(gca, 'fontsize', 12)
            hold on; plot(t_cp_millisecs,exp(CP_OFF),'linewidth',LW);
                ylabel('Gain');
                xlabel('Time [msec]');
                title('Off Parasol Filters');
                plot(zeroline_x, zeroline_y_OFF, 'k');
                ylim([min(.5,min(exp(CP_OFF(:)))), max(1.5,max(exp(CP_OFF(:))))] ); 
            
              %{  
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
             %}
             
             
        end
%%%%%%%%%%%


        orient tall
        eval(sprintf('print -dpdf %s/%s',savedir, filename))
end

    





if strcmp(Basepars.k_filtermode, 'rk2')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    f_opt = opt_param.f;
    %    exitflag = opt_param.exitflag;
    %   figure(1), clf
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    
    
    MU = p_opt(Z.MU);
    K1 = p_opt(Z.SPACE1)*(p_opt(Z.TIME1)');
    K2 = p_opt(Z.SPACE2)*(p_opt(Z.TIME2)');
    K  = K1 + K2;
    PS = ps_basis * p_opt(Z.PS);
    
    %%%% 
    
    clf
    subplot(3,2,[1 2])
    set(gca, 'fontsize', 12)
    c = 0;
    c = c+1;
    text(0, 1-0.1*c,sprintf('%s Fit by %s',info.exp_nm, info.fittype))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Cone Model: %s',info.conemodel))
    c = c+1;
    text(0, 1-0.1*c,sprintf('%s  Cell ID: %4d', info.ctype,info.cid))
    c = c+1;
    text(0, 1-0.1*c,sprintf('Tonic drive with gray stim: %d hz',round(exp(MU)))  )
    c = c+1;
	text(0, 1-0.1*c,sprintf('Optimum fmin: %1.3e',opt_param.f))
	c = c+1;
    text(0, 1-0.1*c,sprintf('fminunc exit flag: %d',opt_param.exitflag))
	axis off
    
    
    
    
    


    %%%%%%%%%%%%%%%%%%%%
    %%%  K stuff %%%%%%%%
    space   = sqrt(length(Z.SPACE1));
    frames  = length(Z.TIME1);
    [~,mx_idx]        = max(abs(K(:)));
	[mx_x,mx_y,mx_fr] = ind2sub([space,space,frames],mx_idx);    
    Kxyt = reshape(K,[space,space,frames]);    
    Kxt  = squeeze(Kxyt(:,mx_y,:));
    Kxy = squeeze(Kxyt(:,:,mx_fr));
    clear mx_x mx_y mx_idx K1 K2
    
	[U,S,V]  = svd(K);
    S        = diag(S);
    spfilter = ( S(1)*U(:,1)*V(5,1) ) / norm( S(1)*U(:,1)*V(5,1) ) ;
    
	subplot(3,2,3)  % label in 50 msec intervals
    set(gca, 'fontsize', 12)
    dur         = 1000*frames*tstim;
    msec_tick   = 50:50:dur;
    frame_tick  = floor(msec_tick / (1000*tstim));
    pixel_tick  = 3:3:space;
	col_ax = [min(K(:)),max(K(:))];
    imagesc(Kxt',col_ax);
    set(gca,'xtick',pixel_tick);
    xlabel('Space [pixel]')
    set(gca, 'ytick', frame_tick, 'yticklabel',msec_tick)
    ylabel('Time [msec]')
    title('Stimulus Filter 1D-Space by Time')  
    colorbar
	subplot(3,2,4)
    set(gca, 'fontsize', 12)
	imagesc(reshape(spfilter,[klen,klen])); colormap jet
	axis image
	set(gca,'xtick',pixel_tick);
	set(gca,'ytick',pixel_tick);
    xlabel('Space [pixel]')
    ylabel('Space [pixel]')
    mx_time = round(1000*tstim*mx_fr);
	title(sprintf('Peak Spatial Filter (%dth frame, ~%d msec)',mx_fr,mx_time))
    colorbar
%%

    LW = 2;
       
    subplot(323)
    set(gca, 'fontsize', 12)
    timebins = length(PS);
    t_ps = 1000*(0:dt:(timebins-1)*dt);
	plot(t_ps,exp(PS),'linewidth',LW)
	hold on
	plot([0,t_ps(end)],[1,1],'k--') % black reference line 
    ylabel('Gain')
    xlabel('Time [msec]')
    title('Post-Spike Filter')
    xlim([0,t_ps(end)])

        

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
           % CP_OFF0 = zeros( max(size(Basepars.cp_basis)) , n_off );
           % CP_ON0  = zeros( max(size(Basepars.cp_basis)) , n_on);

            for i = 1 : n_off
                cid        = Basepars.OFF_Neighbor(i) ;
                n_idx      = find(Basepars.cp_Neighbors == cid);
                paramshift = min(Basepars.paramind.CP) ;

                %%%%%%%%%%%% 
                paramind     = [ (paramshift + (n_idx-1)*Basepars.cp_filternumber) : (paramshift + n_idx*Basepars.cp_filternumber -1) ];
                CP_OFF(:,i)  = Basepars.cp_basis * (p_opt(paramind ));
             %   CP_OFF0(:,i) = Basepars.cp_basis * (p0(paramind));
            end
            for i = 1 : n_on
                cid        = Basepars.ON_Neighbor(i) ;
                n_idx      = find(Basepars.cp_Neighbors == cid);
                paramshift = min(Basepars.paramind.CP) ;

                %%%%%%%%%%%% 
                paramind     = [ (paramshift + (n_idx-1)*Basepars.cp_filternumber) : (paramshift + n_idx*Basepars.cp_filternumber -1) ];
                CP_ON(:,i)  = Basepars.cp_basis * (p_opt(paramind ));
            %    CP_ON0(:,i) = Basepars.cp_basis * (p0(paramind));
            end
            
            
             zeroline_y_ON = linspace( min(exp(CP_ON(:))), max(exp(CP_ON(:))), 100);
             zeroline_y_OFF = linspace( min(exp(CP_ON(:))), max(exp(CP_OFF(:))), 100);
             zeroline_x = zeros(size(zeroline_y_ON));
             LW = 1.3;
             subplot(3,2,5); set(gca, 'fontsize', 12);
             hold on; plot(t_cp_millisecs,exp(CP_ON),'linewidth',LW);   
                ylabel('Gain');
                xlabel('Time [msec]');
                title('On Parasol Filters');
                plot(zeroline_x, zeroline_y_ON, 'k');
                ylim([min(.5,min(exp(CP_ON(:)))), max(1.5,max(exp(CP_ON(:))))] ); 
               
                
            subplot(3,2,6); set(gca, 'fontsize', 12)
            hold on; plot(t_cp_millisecs,exp(CP_OFF),'linewidth',LW);
                ylabel('Gain');
                xlabel('Time [msec]');
                title('Off Parasol Filters');
                plot(zeroline_x, zeroline_y_OFF, 'k');
                ylim([min(.5,min(exp(CP_OFF(:)))), max(1.5,max(exp(CP_OFF(:))))] ); 
            
              %{  
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
             %}
             
             
        end
%%%%%%%%%%%


        orient tall
        eval(sprintf('print -dpdf %s/%s',savedir, filename))
       
end

end
%%%%%
%%


    

%%%%%%%%%
%%%%%% Frame by Frame view %%%%%%%%%%
%{
figure;
detailedfilter_dir = sprintf('%s/FilterDetails',savedir);
if ~isdir(detailedfilter_dir), mkdir(detailedfilter_dir); end
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

