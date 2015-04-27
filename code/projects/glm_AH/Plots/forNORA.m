%%% 2015-02-03
%%% Plot Template for NORA to match AKHeitman plots made for EJs talk

% General Housekeeping and Parameters
clf; hold on; 
set(gca,'fontsize',12);
axis square
MS = 40;
MS_ON = 20;

%%%!!! %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THE LIMITS HERE MAY VARY DEPENDING ON YOUR DATA
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
    if i_exp == 1, colorstring = 'r'; end
    if i_exp == 2, colorstring = 'g'; end
    if i_exp == 3, colorstring = 'b'; end
    
    % Celltype 1 is ON-Parasol,Celltype 2 is OFF-Parasol
    for i_celltype = 1:2
        
        
        %%%!!! %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% NEED TO FILL IN YOUR DATA HERE !!!!! %%%%
        %%%% Vals is a vector of scores, each entry represents a single
        %%%% cell%%%
        if i_celltype == 1, vals_WN = WN_ONP_scores;  vals_NSEM = NSEM_ONP_scores;  end
        if i_celltype == 1, vals_WN = WN_OFFP_scores; vals_NSEM = NSEM_OFFP_scores; end
        
        % Truncate skewed values .. just to make a easier to interpret plot
        vals_WN(find(vals_WN <= minval)) = minval;
        vals_NSEM(find(vals_NSEM <= minval)) = minval;
        
        % Colored Dots
        plotstring = sprintf('%s.',colorstring);
        
        % Plot WN vs. NSEM comparison
        plot(vals_WN, vals_NSEM, plotstring, 'markersize', MS);
        
        % If ON-Parasol indicate with white dot
        if i_celltype == 1, plot(vals_WN, vals_NSEM, 'w.', 'markersize', MS_ON); end
        
    end
end
        
        
        
        
        



