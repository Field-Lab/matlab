function subR_plotfittedNL(fittedGLM, fittedGLM_preNL)
% Cleaned up AKHeitman 2015-06-24

clf;  
printname = sprintf('DiagNLPlot_%s',fittedGLM.cell_savename);
info    = fittedGLM.cellinfo;
GLMType = fittedGLM_preNL.GLMType;

% text box
subplot(3,1,1)
axis off
set(gca, 'fontsize', 12)
obj_NEW = fittedGLM.rawfit.objective_val;
obj_OLD = fittedGLM_preNL.rawfit.objective_val
optNL_describe  = fittedGLM.NL_Output.note_metric;
optNL_string    = fittedGLM.NL_Output.param_string;

c = 0; offset = 0; delta_c = 1.1;
text(-offset, 1-0.1*c,sprintf('%s: %s %d: %s-Fit (red): POSTNL refit with %s',...
    info.exp_nm, info.celltype,info.cid, GLMType.fit_type, fittedGLM.nonlinearity), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('Red is original GLM, Blue includes Postfilter Nonlinearity'))
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('Objval PctChange %d: from %1.2e to %1.2e',...
   round(100*(obj_NEW-obj_OLD)/obj_OLD),obj_OLD,obj_NEW), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('%s', optNL_describe), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('%s', optNL_string), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('Original Fit: %s',GLMType.fitname), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('NL fit Computated at %s',datestr(clock)), 'interpreter','none')
c = c + delta_c; 
text(-offset, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')), 'interpreter','none' );



% plotting non-linearities
x1 = sort(fittedGLM.stimtransform.normalized_filteroutput_fit); 
y1 = sort(fittedGLM.stimtransform.cif_rawGLM_fit);
y2 = sort(fittedGLM.stimtransform.cif_withNL_fit);

LW = 2;
subplot(3,3,4); set(gca, 'fontsize', 10); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([-4,4]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity - Central Portion')

subplot(3,3,5); set(gca, 'fontsize', 10); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([ max(min(x1),-10),0]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity: Inhibitory Portion')

subplot(3,3,6); set(gca, 'fontsize', 10); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([0,min(10,max(x1))]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity: Excitatory Portion')


% plot rasters
subplot(3,1,3); set (gca, 'fontsize',10)
secs     = 8;
dt = fittedGLM.t_bin;
bins     = 120 * 8 * fittedGLM.bins_per_frame;
rec_rast = fittedGLM.xvalperformance.rasters.recorded(:,1:bins);
glm_rast = fittedGLM_preNL.xvalperformance.rasters.glm_sim(:,1:bins);  
NL_rast  = fittedGLM.xvalperformance.rasters.glm_sim(:,1:bins); 
trials   = size(rec_rast,1);
time     = dt*[1:bins];
xlim([0, ceil(time(end))]);
ylim([0 , 3*trials]); hold on
for i_trial = 1:trials
    rec1 = time(find(rec_rast(i_trial,:)));
    glm1 = time(find(glm_rast(i_trial,:)));
    NL1  = time(find(NL_rast(i_trial,:)));
    % Plot the raster
    plot(rec1, i_trial, 'k.')


    yshift = i_trial;
    if length(glm1) < 4*length(rec1) 
        if length(glm1) > 0
            plot(glm1, yshift + trials, 'r.')
        end
    end
    if length(NL1) < 4*length(rec1) 
        if length(NL1) > 0
            plot(NL1, yshift + 2*trials, 'b.')
        end
    end
end
xlabel('seconds'); ylabel('trials')
orient landscape
eval(sprintf('print -dpdf %s.pdf',printname))
end