% AKHEITMAN 2015-07-16
% Define convergence without ignoring improvement in first half
restoredefaultpath; clear
eval(sprintf('load allcells.mat'))
eval(sprintf('load Train_Conv.mat'))


exps = [1] % 2 3 4];
fits = [1 2];

load_scores    = true;
compute_scores = true;
plotting = true;

if load_scores
quartiles = [0 .25 .55 .75 1];
all_percentiles = [0 .5 .15 .25 .35 .45 .55 .65 .75 .85 .95 1]; 
quartile_indices = [3 6 8 11];

aggregated_scores = allcells;
for i_exp = exps

    aggregated_scores{i_exp}.quartiles = quartiles;
    aggregated_scores{i_exp}.all_percentiles = all_percentiles;
    exp_nm = allcells{i_exp}.exp_nm;
    
    % Initialize Structure
    for i_celltype = [1,2]
        if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
        aggregated_scores{i_exp}.celltype{i_celltype}.rawscores_WN   = cell(length(cellgroup),1);
        aggregated_scores{i_exp}.celltype{i_celltype}.rawscores_NSEM = cell(length(cellgroup),1);
    end 
    
    for i_celltype = 1:2
        if i_celltype == 1, ctype = 'ONPar'; cellgroup = allcells{i_exp}.ONP;   end
        if i_celltype == 2, ctype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
        
         for i_fit = fits
            % Pull the scores
            if i_fit == 1
                load_dir   = sprintf('nostim_onlyMU/WN_mapPRJ/%s',exp_nm);
                rawscores = aggregated_scores{i_exp}.celltype{i_celltype}.rawscores_WN; 
            end
            if i_fit == 2
                load_dir   = sprintf('nostim_onlyMU/NSEM_mapPRJ/%s',exp_nm);
                rawscores = aggregated_scores{i_exp}.celltype{i_celltype}.rawscores_NSEM;
            end
            
            if i_fit == 1 && i_celltype == 1
                previous_fits = FullTrain_Conv{i_exp}.test_train_ONP.WN;
            end
            if i_fit == 1 && i_celltype == 2
                previous_fits = FullTrain_Conv{i_exp}.test_train_OFFP.WN;
            end
            if i_fit == 2 && i_celltype == 1
                previous_fits = FullTrain_Conv{i_exp}.test_train_ONP.NSEM;
            end
            if i_fit == 2 && i_celltype == 2
                previous_fits = FullTrain_Conv{i_exp}.test_train_OFFP.NSEM;
            end
            
            for i_cell = 1:length(cellgroup)
                cid = cellgroup(i_cell);
                cell_savename = sprintf('%s_%d',ctype,cid);
                display(sprintf('working on Piece: %s, Cell: %s', exp_nm, cell_savename));
                eval(sprintf('load %s/%s.mat fittedGLM', load_dir, cell_savename));
                new_pt        = fittedGLM.rawfit.objective_val;
                raw_quartiles = [new_pt; [previous_fits{i_cell}.objval(quartile_indices)]];
                raw_allpts    = [new_pt; [previous_fits{i_cell}.objval]];
                rawscores{i_cell}.raw_quartiles = raw_quartiles;
                rawscores{i_cell}.raw_allpts    = raw_allpts;
                clear fittedGLM
            end
            
            if i_fit == 1
                aggregated_scores{i_exp}.celltype{i_celltype}.rawscores_WN = rawscores;
            end
            if i_fit == 2
                aggregated_scores{i_exp}.celltype{i_celltype}.rawscores_NSEM = rawscores;
            end
            
         end
    end    
    
end
eval(sprintf('save convergence_aggscores.mat aggregated_scores'));

end
%%
if compute_scores  
    eval(sprintf('load convergence_aggscores.mat'));
for i_exp = 1:4
    
    
    xvals = aggregated_scores{i_exp}.quartiles;
    decay_constants = 2.^[-1:.1:5];
    aggregated_scores{i_exp}.decay_constants = decay_constants;
    % Initialize Structure
    for i_celltype = 1:2
        if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
        aggregated_scores{i_exp}.celltype{i_celltype}.fittedexp_WN   = NaN(length(cellgroup),1);
        aggregated_scores{i_exp}.celltype{i_celltype}.fittedexp_NSEM = NaN(length(cellgroup),1);
        
        aggregated_scores{i_exp}.celltype{i_celltype}.halfratio_WN   = NaN(length(cellgroup),1);
        aggregated_scores{i_exp}.celltype{i_celltype}.halfratio_NSEM = NaN(length(cellgroup),1);
    end 
    
    
    
    for i_celltype = 1:2
        if i_celltype == 1, celltype = 'ONPar';  cellgroup = allcells{i_exp}.ONP;  end
        if i_celltype == 2, celltype = 'OFFPar'; cellgroup = allcells{i_exp}.OFFP; end
         for i_fit = 1:2
            % Pull the scores
            if i_fit == 1
                fittedexp = aggregated_scores{i_exp}.celltype{i_celltype}.fittedexp_WN; 
                halfratio = aggregated_scores{i_exp}.celltype{i_celltype}.halfratio_WN;
                rawscores = aggregated_scores{i_exp}.celltype{i_celltype}.rawscores_WN; 
            end
            if i_fit == 2
                fittedexp = aggregated_scores{i_exp}.celltype{i_celltype}.fittedexp_NSEM; 
                rawscores = aggregated_scores{i_exp}.celltype{i_celltype}.rawscores_NSEM;
                halfratio = aggregated_scores{i_exp}.celltype{i_celltype}.halfratio_NSEM;
            end

            for i_cell = 1:length(cellgroup)
                rawvals = rawscores{i_cell}.raw_quartiles;
                rawvals = (rawvals - rawvals(1));
                rawvals = rawvals / abs(rawvals(1)-rawvals(end))  + 1;
                
                halfratio(i_cell) = rawvals(3);
                
                % ASSUME EXPONENTIAL(-1/tau)  FIND CORRESPONDING TOP POINT
                error_vec = zeros(1,length(decay_constants));
                for i_decay   = 1:length(decay_constants)
                    net_decay = decay_constants(i_decay);
                    final_val = exp(-net_decay);
                    
                    yvals     = (1-final_val)*rawvals + final_val;
                    yvals_fit = exp(-net_decay.*(xvals));
                    error_vec(i_decay) = norm(yvals_fit' - yvals);
                end
                [error, index] = min(error_vec);
                fittedexp(i_cell) = decay_constants(index);
            end
            if i_fit == 1
                aggregated_scores{i_exp}.celltype{i_celltype}.fittedexp_WN   = fittedexp;
                aggregated_scores{i_exp}.celltype{i_celltype}.halfratio_WN   = halfratio;
            end
            if i_fit == 2
                aggregated_scores{i_exp}.celltype{i_celltype}.fittedexp_NSEM = fittedexp;
                aggregated_scores{i_exp}.celltype{i_celltype}.halfratio_NSEM = halfratio;
            end
         end
    end
    
    
end

eval(sprintf('save fitted_aggscores.mat aggregated_scores'));
end    


%%



if plotting
 eval(sprintf('load fitted_aggscores.mat'));
    
    
for i_type = 1:2    
   if i_type == 1 
       type = 'fittedexp';  
       low_lim  = 0;
       high_lim = 20;
   end
   if i_type == 2
        type    = 'halfratio';    
        low_lim = 0;
        high_lim = .1;
   end
   colors = {'r.','g.','b.','c.'};
    for i_exp = 1:4
        
        clf


        colorstring = colors{i_exp};
        for i_celltype = [1:2]
            subplot(1,2,i_celltype)

            MS_A = 12;
            MS_B = 30; 
            hold on; set(gca, 'fontsize', 12)
            axis square
            xlim([low_lim, high_lim]);
            ylim([low_lim, high_lim]);
            plot(linspace(low_lim, high_lim,100),linspace(low_lim, high_lim,100),'k');



            if strcmp(type, 'halfratio')
                scores1 = aggregated_scores{i_exp}.celltype{i_celltype}.halfratio_WN;
                scores2 = aggregated_scores{i_exp}.celltype{i_celltype}.halfratio_NSEM;
                title('Improvement 2nd half / Improvement 1st half');
            elseif strcmp(type, 'fittedexp')
                scores1 = aggregated_scores{i_exp}.celltype{i_celltype}.fittedexp_WN;
                scores2 = aggregated_scores{i_exp}.celltype{i_celltype}.fittedexp_NSEM;
                title('Net Decay Time Constants');
            end

            scores1(find(scores1<=low_lim)) = low_lim;
            scores2(find(scores2<=low_lim)) = low_lim;
            scores1(find(scores1>=high_lim)) = high_lim;
            scores2(find(scores2>=high_lim)) = high_lim;

            plot(scores1, scores2, colorstring, 'markersize', MS_B);
            if i_celltype == 1
                plot(scores1, scores2, 'w.', 'markersize', MS_A);
            elseif i_celltype == 2
                plot(scores1, scores2, 'k.', 'markersize', MS_A);
            end 
            
            xlabel('White Noise')
            ylabel('NSEM')
            
        end    
        orient landscape
        if strcmp(type, 'halfratio')
            eval(sprintf('print -dpdf halfratio_%d.pdf',i_exp))
        elseif strcmp(type, 'fittedexp')
            eval(sprintf('print -dpdf decaytime_%d.pdf',i_exp))
        end
    end
end


end
 
        
       
        
        
        


