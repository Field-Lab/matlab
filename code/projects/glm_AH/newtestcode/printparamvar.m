function printparamvar(fittedGLM,manual_search,printname)

info    = fittedGLM.cellinfo;
GLMType = fittedGLM.GLMType;
dt = fittedGLM.t_bin;
tstim = fittedGLM.bins_per_frame * dt;
performance = fittedGLM.xvalperformance;
if GLMType.TonicDrive
    MU = fittedGLM.linearfilters.TonicDrive.Filter;
end
if GLMType.PostSpikeFilter
    PS = fittedGLM.linearfilters.PostSpike.Filter;
end

K        = fittedGLM.linearfilters.Stimulus.Filter;
K_time1  = fittedGLM.linearfilters.Stimulus.time_rk1; 
K_space1 = fittedGLM.linearfilters.Stimulus.space_rk1; 
    
clf
subplot(5,1,1)
axis off
set(gca, 'fontsize', 12)
c = 0;
text(-.1, 1-0.1*c,sprintf('%s: %s %d: %s-Fit',info.exp_nm, info.celltype,info.cid, GLMType.fit_type))
c = c + 1.5;
text(-.1, 1-0.1*c,sprintf('Xval performance, normed by unconditioned optimal performance (UOP) %1.3e',performance.glm_normedbits))
 c = c + 2;
text(-.1, 1-0.1*c,sprintf('Fit Type: %s',GLMType.fitname), 'interpreter','none')
c = c + 1.5;
if GLMType.TonicDrive
    text(-.1, 1-0.1*c,sprintf('Tonic drive with gray stim: %1.2e hz',round(exp(MU)))  )
    c = c + 1.5;
end

text(-.1, 1-0.1*c,sprintf('Optimum fmin: %1.4e',fittedGLM.rawfit.objective_val))
c = c + 1.5;
text(-.1, 1-0.1*c,sprintf('Raw xval bitsperspike %1.3e',performance.logprob_glm_bpspike))
c = c + 1.2;
text(-.1, 1-0.1*c,sprintf('Computated at %s',fittedGLM.fit_time));





nonlinpar = manual_search.nonlinpar;
lin_index = manual_search.lin_index;

obj_val = - (manual_search.objective_value - manual_search.objective_value(lin_index) );
obj_val = obj_val / (abs(manual_search.objective_value(lin_index)));

maxval = max(obj_val);
max_index = find(obj_val == maxval);

if length(max_index > 1)
    max_index = max_index(1);
end

if maxval == 0
    maxval = mean(abs(obj_val));
    max_index = lin_index;
end

plotobj_val = obj_val;
for i_ind = 1:length(obj_val(:));
    if obj_val(i_ind) < -maxval;
        plotobj_val(i_ind) = -maxval;
    end
end

if manual_search.dimension == 1
    subplot(5,1,[2,3,4,5])
    plot(linspace(nonlinpar(1),nonlinpar(end),100), zeros(1,100),'k'); hold on;
    plot(nonlinpar , plotobj_val, 'r', 'linewidth', 2);
    plot(nonlinpar(max_index), maxval, 'k*', 'markersize', 20); 
    ylim([-1*maxval, maxval]); xlim([nonlinpar(1), nonlinpar(end)]);
    
    ylabel('Improvement from Linear / Linear Score');
    xlabel('Non-linearity param');
end

if manual_search.dimension == 2
    subplot(5,1,[2,3,4,5])
    if strcmp(manual_search.NL, 'piece_linear_shiftmean')
        val_1 = manual_search.nonlinpar_ratio; y_label = 'Increment to Decrement ratio';
        val_2 = manual_search.nonlinpar_shift; x_label = 'Shift the Zero Point';
        header = 'Moving Hinge';
    end
    linpar = nonlinpar(lin_index,:);
    lin_1  = find(val_1 == linpar(1));
    lin_2  = find(val_2 == linpar(2));
    
    maxpar = nonlinpar(max_index, :);
    max_1  = find(val_1 == maxpar(1));
    max_2  = find(val_2 == maxpar(2));
    
    b = reshape(plotobj_val, [length(val_1),length(val_2)]);  
    set(gca, 'fontsize', 12); 
    imagesc(b); colorbar; colorbar('fontsize',12); colormap gray
    hold on;
    
    %%%  Just an annoying flip due to the way I'm plotting %%%
    % y_tick goes to value 1 %
    % x_tick goes to value 2 %
    y_tick  = round(length(val_1) * [.2 .4 .6 .8 1]);
    y_ticklabel = val_1(y_tick);
    
    
    x_tick  = round(length(val_2) * [.2 .4 .6 .8 1]);
    x_ticklabel = val_2(x_tick);
    
    set(gca, 'xtick', x_tick, 'xticklabel',x_ticklabel)
    set(gca, 'ytick', y_tick, 'yticklabel',y_ticklabel)
    xlabel(x_label); ylabel(y_label);
    
    plot(lin_2, lin_1, 'r.','markersize',30);
    plot(max_2, max_1, 'g.','markersize',30);
    
    title(sprintf('%s: Linear (Red Point), Optimum (Green)', header));
end


if manual_search.dimension > 2
    
    x_dim = linspace(0,1,1000);
    if strcmp(manual_search.NL, 'polynomial_order3_part5')
        opt_param = fittedGLM.pt_nonlinearity_param;  
    end

    [ignore,list] = sort(obj_val, 'descend');
    
    topfive    =list(1:5);
    bottomfive =list((end-4):end); 
    
    top_y    = zeros(5,length(x_dim));
    bottom_y = zeros(5,length(x_dim));
    
    tick_marks = [0, .25, .5, .75, 1];
    
    
    for i_count = 1:5
        NL_par_top    = nonlinpar(topfive(i_count),:);
        NL_par_bot    = nonlinpar(bottomfive(i_count),:);
        if strcmp(manual_search.NL, 'polynomial_order3_part5')
            y1 = NL_par_top(1)*(x_dim.^(1/1)) + ...
                 NL_par_top(2)*(x_dim.^(2/1)) +  NL_par_top(3)*(x_dim.^(1/2)) + ...
                 NL_par_top(4)*(x_dim.^(3/1)) +  NL_par_top(5)*(x_dim.^(1/3)) ;
            y2 = NL_par_bot(1)*(x_dim.^(1/1)) + ...
                 NL_par_bot(2)*(x_dim.^(2/1)) +  NL_par_bot(3)*(x_dim.^(1/2)) + ...
                 NL_par_bot(4)*(x_dim.^(3/1)) +  NL_par_bot(5)*(x_dim.^(1/3)) ;
        end
        if strcmp(manual_search.NL, 'polynomial_order5_part4')
            y1 = NL_par_top(1)*(x_dim.^(1/1)) + ...
                 NL_par_top(2)*(x_dim.^(2/1)) +  NL_par_top(3)*(x_dim.^(1/2)) + ...
                 NL_par_top(4)*(x_dim.^(3/1)) +  NL_par_top(5)*(x_dim.^(1/3))  + ...
                 NL_par_top(6)*(x_dim.^(4/1)) +  NL_par_top(7)*(x_dim.^(1/4))  + ...
                 NL_par_top(8)*(x_dim.^(5/1)) +  NL_par_top(9)*(x_dim.^(1/5)) ;
                 
            y2 = NL_par_bot(1)*(x_dim.^(1/1)) + ...
                 NL_par_bot(2)*(x_dim.^(2/1)) +  NL_par_bot(3)*(x_dim.^(1/2)) + ...
                 NL_par_bot(4)*(x_dim.^(3/1)) +  NL_par_bot(5)*(x_dim.^(1/3))  + ...
                 NL_par_bot(6)*(x_dim.^(4/1)) +  NL_par_bot(7)*(x_dim.^(1/4))  + ...
                 NL_par_bot(8)*(x_dim.^(5/1)) +  NL_par_bot(9)*(x_dim.^(1/5))  ;
        end
        
        if strcmp(manual_search.NL, 'piecelinear_fourpiece_eightlevels')
            
            for i_type = 1:2
                if i_type == 1, nLpar = NL_par_top; end
                if i_type == 2, nLpar = NL_par_bot; end
                    
                slope1  = nLpar(1);    
                slope2  = nLpar(2);    
                slope3  = nLpar(3);
                slope4  = nLpar(4);
                offset1 = 0;
                offset2 = .25* slope1;
                offset3 = .25*(slope1+slope2);
                offset4 = .25*(slope1+slope2+slope3);
                ydummy  = x_dim;
                ydummy(1:250) = slope1 * (ydummy(  1:250)  - .00) + offset1;
                ydummy(251:500) = slope2 * (ydummy(251:500)  - .25) + offset2;
                ydummy(501:750) = slope3 * (ydummy(501:750)  - .50) + offset3;
                ydummy(751:1000) = slope4 * (ydummy(751:1000) - .75) + offset4;
                
                if i_type == 1, y1 = ydummy; end
                if i_type == 2, y2 = ydummy; end
                
                    
            end
        end
        
        
        top_y(i_count,:)    = y1;
        bottom_y(i_count,:) = y2;
    end
    
    subplot(5,2,[3,5,7,9])
    set(gca,'fontsize', 10);
    set(gca,'xtick',tick_marks);
    set(gca,'ytick',tick_marks);
    xlim([0 1]);
    ylim([0,1]);
    xlabel('Stimulus Values before Non-Lin');
    ylabel('Non-Lin Transformed Values');
    title('Magenta = Best, Red = next 5, Blue = last 5');
    hold on
    for i_count = 1:5
        plot(x_dim, top_y(i_count,:), 'r');
        plot(x_dim, bottom_y(i_count,:),'b');
    end
    plot(x_dim, top_y(1,:) , 'm', 'linewidth', 3);
    
        
    
    
    quartile        = round(length(list)/4);
    topquartile     = list(1:quartile);
    bottomquartile  = list((end-quartile+1):end);
    middlequartiles = list(quartile:(end-quartile));
    
    topq_y    = zeros(quartile, length(x_dim));
    bottomq_y = zeros(quartile, length(x_dim));
    for i_count = 1:quartile
        NL_par_top    = nonlinpar(topquartile(i_count),:);
        NL_par_bot    = nonlinpar(bottomquartile(i_count),:);
        if strcmp(manual_search.NL, 'polynomial_order3_part5')
            y1 = NL_par_top(1)*(x_dim.^(1/1)) + ...
                 NL_par_top(2)*(x_dim.^(2/1)) +  NL_par_top(3)*(x_dim.^(1/2)) + ...
                 NL_par_top(4)*(x_dim.^(3/1)) +  NL_par_top(5)*(x_dim.^(1/3)) ;
            y2 = NL_par_bot(1)*(x_dim.^(1/1)) + ...
                 NL_par_bot(2)*(x_dim.^(2/1)) +  NL_par_bot(3)*(x_dim.^(1/2)) + ...
                 NL_par_bot(4)*(x_dim.^(3/1)) +  NL_par_bot(5)*(x_dim.^(1/3)) ;    
        end
        
        if strcmp(manual_search.NL, 'polynomial_order5_part4')
            y1 = NL_par_top(1)*(x_dim.^(1/1)) + ...
                 NL_par_top(2)*(x_dim.^(2/1)) +  NL_par_top(3)*(x_dim.^(1/2)) + ...
                 NL_par_top(4)*(x_dim.^(3/1)) +  NL_par_top(5)*(x_dim.^(1/3))  + ...
                 NL_par_top(6)*(x_dim.^(4/1)) +  NL_par_top(7)*(x_dim.^(1/4))  + ...
                 NL_par_top(8)*(x_dim.^(5/1)) +  NL_par_top(9)*(x_dim.^(1/5)) ;
                 
            y2 = NL_par_bot(1)*(x_dim.^(1/1)) + ...
                 NL_par_bot(2)*(x_dim.^(2/1)) +  NL_par_bot(3)*(x_dim.^(1/2)) + ...
                 NL_par_bot(4)*(x_dim.^(3/1)) +  NL_par_bot(5)*(x_dim.^(1/3))  + ...
                 NL_par_bot(6)*(x_dim.^(4/1)) +  NL_par_bot(7)*(x_dim.^(1/4))  + ...
                 NL_par_bot(8)*(x_dim.^(5/1)) +  NL_par_bot(9)*(x_dim.^(1/5))  ;
        end
        if strcmp(manual_search.NL, 'piecelinear_fourpiece_eightlevels')
        
            for i_type = 1:2
                if i_type == 1, nLpar = NL_par_top; end
                if i_type == 2, nLpar = NL_par_bot; end
                    
                slope1  = nLpar(1);    
                slope2  = nLpar(2);    
                slope3  = nLpar(3);
                slope4  = nLpar(4);
                offset1 = 0;
                offset2 = .25* slope1;
                offset3 = .25*(slope1+slope2);
                offset4 = .25*(slope1+slope2+slope3);
                ydummy  = x_dim;
                ydummy(1:250) = slope1 * (ydummy(  1:250)  - .00) + offset1;
                ydummy(251:500) = slope2 * (ydummy(251:500)  - .25) + offset2;
                ydummy(501:750) = slope3 * (ydummy(501:750)  - .50) + offset3;
                ydummy(751:1000) = slope4 * (ydummy(751:1000) - .75) + offset4;
                
                if i_type == 1, y1 = ydummy; end
                if i_type == 2, y2 = ydummy; end
                
                    
            end
        end
        
        
        topq_y(i_count,:)    = y1;
        bottomq_y(i_count,:) = y2;
    end
        
    middle_y = zeros(size(middlequartiles,1) , length(x_dim));
    for i_count = 1:size(middlequartiles,1)
        NL_par_mid    = nonlinpar(middlequartiles(i_count),:);
        if strcmp(manual_search.NL, 'polynomial_order3_part5')
            y1 = NL_par_mid(1)*(x_dim.^(1/1)) + ...
                 NL_par_mid(2)*(x_dim.^(2/1)) +  NL_par_mid(3)*(x_dim.^(1/2)) + ...
                 NL_par_mid(4)*(x_dim.^(3/1)) +  NL_par_mid(5)*(x_dim.^(1/3)) ;
        end
        if strcmp(manual_search.NL, 'polynomial_order5_part4')
            y1 = NL_par_mid(1)*(x_dim.^(1/1)) + ...
                 NL_par_mid(2)*(x_dim.^(2/1)) +  NL_par_mid(3)*(x_dim.^(1/2)) + ...
                 NL_par_mid(4)*(x_dim.^(3/1)) +  NL_par_mid(5)*(x_dim.^(1/3))  + ...
                 NL_par_mid(6)*(x_dim.^(4/1)) +  NL_par_mid(7)*(x_dim.^(1/4))  + ...
                 NL_par_mid(8)*(x_dim.^(5/1)) +  NL_par_mid(9)*(x_dim.^(1/5)) ;
                 
            
        end
        if strcmp(manual_search.NL, 'piecelinear_fourpiece_eightlevels')
        
            for i_type = 1
                if i_type == 1, nLpar = NL_par_mid; end
                    
                slope1  = nLpar(1);    
                slope2  = nLpar(2);    
                slope3  = nLpar(3);
                slope4  = nLpar(4);
                offset1 = 0;
                offset2 = .25* slope1;
                offset3 = .25*(slope1+slope2);
                offset4 = .25*(slope1+slope2+slope3);
                ydummy  = x_dim;
                ydummy(1:250) = slope1 * (ydummy(  1:250)  - .00) + offset1;
                ydummy(251:500) = slope2 * (ydummy(251:500)  - .25) + offset2;
                ydummy(501:750) = slope3 * (ydummy(501:750)  - .50) + offset3;
                ydummy(751:1000) = slope4 * (ydummy(751:1000) - .75) + offset4;
                
                if i_type == 1, y1 = ydummy; end
                
                    
            end
        end
        
        
        middle_y(i_count,:)    = y1;        
    end
    
    subplot(5,2,[4,6,8,10])
    set(gca,'fontsize', 10);
    set(gca,'xtick',tick_marks);
    set(gca,'ytick',tick_marks);
    xlim([0 1]);
    ylim([0,1]);
    xlabel('Stimulus Values before Non-Lin');
    ylabel('Non-Lin Transformed Values');
    title('Red = Top Quartile, Black = Middle, Blue = Bottom Quartile');
    hold on
    for i_count = 1:size(middlequartiles,1)
        plot(x_dim, middle_y(i_count,:), 'k');
    end
    for i_count = 1:quartile
        plot(x_dim, topq_y(i_count,:), 'r');
        plot(x_dim, bottomq_y(i_count,:),'b');
    end

    plot(x_dim, top_y(1,:) , 'm', 'linewidth', 6);
    
   
end

    
    
    

orient landscape
eval(sprintf('print -dpdf %s.pdf',printname))
%eval(sprintf('print -dpng %s.png',printname))
end
    
