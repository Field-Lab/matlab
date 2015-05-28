%%% 2015-02-03
%%% Plot Template for NORA to match AKHeitman plots made for EJs talk
figure; hold on

% Data location(added by Nora)
fit_type='fixedSP_rk1_linear_MU_PS_CP_p8IDp8';
fit_directory=['/Volumes/Analysis/nora/NSEM/GLM_Output/' fit_type '/standardparams'];

% General Housekeeping and Parameters
clf; hold on; 
set(gca,'fontsize',12);
axis square
MS = 40;
MS_ON = 20;

%%%!!! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THE LIMITS HERE MAY VARY DEPENDING ON YOUR DATA
minval = 0; 
maxval = 0.75;
% maxval = .25;
% set(gca,'xtick', [-.75,-.5,-.25,0,.25])
% set(gca,'ytick', [-.75,-.5,-.25,0,.25])

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
    
    % Added by nora to collect data
    cells_NSEM=dir([fit_directory '/NSEM_mapPRJ/' exp '/']);
    cells_WN=dir([fit_directory '/WN_mapPRJ/' exp '/']);
    ON_cell_count=1;
    OFF_cell_count=1;
    WN_ONP_scores=zeros(100,1);
    WN_OFFP_scores=zeros(100,1);
    NSEM_ONP_scores=zeros(100,1);
    NSEM_OFFP_scores=zeros(100,1);
        
    for i_cell=1:length(cells_NSEM)
        if cells_NSEM(i_cell).name(1) == 'O'
            if strcmp(cells_NSEM(i_cell).name(2),'N')
                try
                    load([fit_directory '/WN_mapPRJ/' exp '/' cells_NSEM(i_cell).name])
                    WN_ONP_scores(ON_cell_count)=fittedGLM.xvalperformance.glm_normedbits;
                    load([fit_directory '/NSEM_mapPRJ/' exp '/' cells_NSEM(i_cell).name])
                    NSEM_ONP_scores(ON_cell_count)=fittedGLM.xvalperformance.glm_normedbits;
                    ON_cell_count=ON_cell_count+1;
                catch
                    warn([cells_NSEM(i_cell).name ' does not exist in white noise']);
                end
            else
                try
                    load([fit_directory '/WN_mapPRJ/' exp '/' cells_NSEM(i_cell).name])
                    WN_OFFP_scores(OFF_cell_count)=fittedGLM.xvalperformance.glm_normedbits;
                    load([fit_directory '/NSEM_mapPRJ/' exp '/' cells_NSEM(i_cell).name])
                    NSEM_OFFP_scores(OFF_cell_count)=fittedGLM.xvalperformance.glm_normedbits;
                    OFF_cell_count=OFF_cell_count+1;
                catch
                    warn([cells_NSEM(i_cell).name ' does not exist in white noise']);
                end
            end
        end
    end
    
    % Celltype 1 is ON-Parasol,Celltype 2 is OFF-Parasol
    for i_celltype = 1:2
        
        %%%!!! %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% NEED TO FILL IN YOUR DATA HERE !!!!! %%%%
        %%%% Vals is a vector of scores, each entry represents a single
        %%%% cell%%%
        if i_celltype == 1, vals_WN = WN_ONP_scores(1:ON_cell_count-1);  vals_NSEM = NSEM_ONP_scores(1:ON_cell_count-1);  end
        if i_celltype == 2, vals_WN = WN_OFFP_scores(1:OFF_cell_count-1); vals_NSEM = NSEM_OFFP_scores(1:OFF_cell_count-1); end
        
        % Truncate skewed values .. just to make a easier to interpret plot
        vals_WN(find(vals_WN <= minval)) = minval;
        vals_NSEM(find(vals_NSEM <= minval)) = minval;
        
        % Colored Dots
        plotstring = sprintf('%s.',colorstring);
        
        % Plot WN vs. NSEM comparison
        plot(vals_WN, vals_NSEM, plotstring, 'markersize', MS);
        
        % If ON-Parasol indicate with white dot
        if i_celltype == 1, plot(vals_WN, vals_NSEM, 'w.', 'markersize', MS_ON);end
       
    end
end
        
        
        
        
        



