spec.type = 'REVERSING-SINUSOID';
spec.rgb = {[0.48]  [0.48]  [0.48]};
spec.orientation = 45;
spec.frames = 480;
spec.x_start = 0;
spec.x_end = 800;
spec.y_start = 0;
spec.y_end = 600;
spec.spatial_period = 80;
spec.temporal_period = 30;
spec.spatial_phase = 3;

colormap(gray);
switch spec.type
    case 'REVERSING-SINUSOID'
    [maxframe, tscale] = calc_reversing_grating_frame_intensities(spec);
    for i=1:length(tscale)
        frame_channel1 = maxframe.*tscale(i).*spec.rgb{1}+0.5;
        imagesc(frame_channel1,[0 1]); axis equal off;
        pause(0.01);
    end
    
    case 'DRIFTING-SINUSOID'
    [framesin, framecos, tdrift] = calc_drifting_grating_frame_intensities_jf(spec);
    for i=1:length(tdrift)
        frameshifted = framesin.*cos(tdrift(i)) + framecos.*sin(tdrift(i));
        framescaled = (0.5*frameshifted);
        frame_channel1 = framescaled.*spec.rgb{1} + 0.5;
        imagesc(frame_channel1,[0 1]); axis equal off;
        pause(0.01);
    end
end
    
    
