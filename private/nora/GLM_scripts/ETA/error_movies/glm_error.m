
function [glm_error] = glm_error(fittedGLM)
%%

% Calculate the PSTH
            convolve=100;
            rec_rast = fittedGLM.xvalperformance.rasters.recorded;
            sim_rast = fittedGLM.xvalperformance.rasters.glm_sim;
            trials = size(rec_rast,1);
            for i=1:trials
                PSTH_rec(i,:)=conv(rec_rast(i,:),ones(convolve,1),'same');
                PSTH_sim(i,:)=conv(sim_rast(i,:),ones(convolve,1),'same');
            end
            PSTH = zeros(2, size(PSTH_rec,2));
            PSTH(1,:) = sum(PSTH_rec);
            PSTH(2,:) = sum(PSTH_sim);
            disp('PSTH calculated');
            glm_error = PSTH(1,:) - PSTH(2,:);

end