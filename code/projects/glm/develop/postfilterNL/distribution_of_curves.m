clear; close all; clc
K_vec  = [1/8 .25 .5 1 2 4 8];
X_vec_ofK = [-4 -2 -1 -.5 0 .5 1 2 4];

SP0 = linspace(-4,4,100);
for i_K = 1:length(K_vec)
    figure; LW = 4
    for i_X0 = 1:length(X_vec_ofK)
        K = K_vec(i_K);
        X0 = X_vec_ofK(i_X0)/ K;
        logistic_vals = 1./ (1 + exp(-K * ( SP0 - X0) ) ) + (  1- 1/(1+exp(K*X0)) );
        % plot(sort(SP0), sort(logistic_vals))
        ymax = 1.2* (max(logistic_vals)-1) + 1;
        
        subplot(3,3, i_X0); hold on
        plot(SP0, logistic_vals, 'b','linewidth', LW)
        ylim([0,ymax]);  set(gca, 'ytick', [1]);
        set(gca,'fontsize',10)
        title(sprintf('Offset: %d, Slope: %1.1e', X_vec_ofK(i_X0),K));
        
        
    end
    orient landscape
    eval(sprintf('print -dpdf Logistics_Slope_%1.1e_withoffsets.pdf', K))
end