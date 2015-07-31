
function [glm_error] = glm_error(fittedGLM)
%%
% Calculate the PSTH
            convolve=100;
            rec_rast = fittedGLM.xvalperformance.rasters.recorded;
            sim_rast = fittedGLM.xvalperformance.rasters.glm_sim;
            PSTH = conv2([sum(rec_rast); sum(sim_rast)], gausswin(convolve)', 'same');
            glm_error = PSTH(1,:)./PSTH(2,:);
end