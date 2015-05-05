STA = zeros(40,80,30);
fitframes = 12000; 

for i = 1:length(fitspikes)
    sp_frame = floor(fitspikes(i) * 120);
    if sp_frame > 29 && sp_frame<fitframes 
        STA = STA+fitmovie(:,:,(sp_frame-29):sp_frame);
    end
end