function fittedGLM = save_glm(fittedGLM, testspikes, testmovie, savename)
    fittedGLM.xvalperformance = glm_predict(fittedGLM,testmovie, 'testspikes', testspikes);
    temp = corrcoef(conv(sum(fittedGLM.xvalperformance.rasters.glm_sim), gausswin(100)),conv(sum(fittedGLM.xvalperformance.rasters.recorded), gausswin(100)));
    fittedGLM.xvalperformance.corr = temp(2,1);
    close all
    plotfilters(fittedGLM)
    exportfig(gcf, [savename '_filters.eps'], 'Bounds', 'loose', 'Color', 'rgb');
    close all
    plotrasters(fittedGLM.xvalperformance, fittedGLM)
    exportfig(gcf, [savename '_rasters.eps'], 'Bounds', 'loose', 'Color', 'rgb');
    close all
    save([savename '.mat'], 'fittedGLM');
end
