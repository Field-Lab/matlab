%%% 2015-02-03
%%% Plot Template for NORA to match AKHeitman plots made for EJs talk

% General Housekeeping and Parameters
clf; hold on; 
set(gca,'fontsize',12);
axis square
MS = 40;
MS_ON = 20;

% Data location(added by Nora)
fit_type='fixedSP_rk1_linear_MU_PS_CP_p8IDp8';
fit_directory=['/Volumes/Analysis/nora/NSEM/GLM_Output/' fit_type '/standardparams'];


%%% !!!! EDIT !!! %%%
% THE LIMITS HERE MAY VARY DEPENDING ON YOUR DATA
minval = -.75; 
maxval = .25;
set(gca,'xtick', [-.75,-.5,-.25,0,.25])
set(gca,'ytick', [-.75,-.5,-.25,0,.25])




% Set up plot and unity line
xlim([minval,maxval]);
ylim([minval,maxval]);
unity_line = linspace(minval, maxval,100);
plot(unity_line,unity_line,'k')



% 2012-08-09-3 as experiment 1
% 2012-09-27-3 as experiment 2
% 2013-08-19-6 as experiment 3
% 2013-10-10-0  as experiment 4 (temporarily droppped)
for i_exp = 1:3
    if i_exp == 1, colorstring = 'r'; exp='2012-08-09-3'; end
    if i_exp == 2, colorstring = 'g'; exp='2012-09-27-3'; end
    if i_exp == 3, colorstring = 'b'; exp='2013-08-19-6'; end
    
    cells_NSEM=dir([fit_directory '/NSEM_mapPRJ/' exp '/']);
    cells_WN=dir([fit_directory '/WN_mapPRJ/' exp '/']);
    
    % Celltype 1 is ON-Parasol,Celltype 2 is OFF-Parasol
    for i_celltype = 1:2
        
        %%% !!!! EDIT !!! %%%
        %Need to load cell ids here ONPars or OFFPARS
        cells = cells %
        
        for cid = cells
            
            %%% !!!! EDIT !!! %%%
            %NEED TO LOAD UP Raw Bits Per Spike of glm %
            rawscore_WN   = xval_WN.logprob_glm_bpspike;
            rawscore_NSEM =  xval_NSEM.logprob_glm_bpspike;
            
            % Identify a numbering index in the new raster_scores structure
            if i_celltype == 1, rast_normindex = find(raster_scores{i_exp}.ONP == cid; end
            if i_celltype == 2, rast_normindex = find(raster_scores{i_exp}.OFFP == cid; end
            
            % Updated Bits Per Spike for the Raster predicting itself
            updated_uop_WN   = raster_scores{i_exp}.stim_type{1}.celltype{i_celltype}.scores.values(rast_normindex);
            updated_uop_NSEM = raster_scores{i_exp}.stim_type{2}.celltype{i_celltype}.scores.values(rast_normindex);
            
            % Subtractive Normalization
            score_WN   =    rawscore_WN   - updated_uop_WN; 
            score_NSEM =    rawscore_NSEM - updated_uop_NSEM;
            
            % Truncate if necessary to keep the plot pretty
            if score_WN <= minval, score_WN = minval; end
            if score_NSEM <= minval, socre_NSEM = minval; end
        
            % Colored Dots
            plotstring = sprintf('%s.',colorstring);
        
            % Plot WN vs. NSEM comparison
            plot(score_WN, score_NSEM, plotstring, 'markersize', MS);
        
            % If ON-Parasol indicate with white dot
            if i_celltype == 1, plot(score_WN, score_NSEM, 'w.', 'markersize', MS_ON); end
        end
        
    end
end
        
        
        
        
        



