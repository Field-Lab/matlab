%cid = [ 32 768 2778 4354 5866 7036];  ON 2013-10-10-0


cidOFF  = [ 1 31 737 1328 1341 2959 5447] %2013-08-19-6
for i_cell = 1:length(cidOFF)
    eval(sprintf('load OFFPar_%d.mat', cidOFF(i_cell) ) );
    figure
    plot(fittedGLM.linearfilters_linstim.Stimulus.time_rk1)
    hold on
    plot(fittedGLM.linearfilters_nonlinstim.Stimulus.time_rk1,'r')
    
    xlabel(sprintf('nonlinparam %d', fittedGLM.pt_nonlinearity_param));
end



cidON  = [ 2824 3167 3996 5660 6799]  % 2013-08-19-6
for i_cell = 1:length(cidON)
    eval(sprintf('load ONPar_%d.mat', cidON(i_cell) ) );
    figure
    plot(fittedGLM.linearfilters_linstim.Stimulus.time_rk1)
    hold on
    plot(fittedGLM.linearfilters_nonlinstim.Stimulus.time_rk1,'r')
    
    xlabel(sprintf('nonlinparam %d', fittedGLM.pt_nonlinearity_param));
end