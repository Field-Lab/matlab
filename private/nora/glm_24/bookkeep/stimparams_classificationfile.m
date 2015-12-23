% AKHeitman 2014-04-02
% Dumb way to read in Noise stimuli parameters

function params = stimparams_classificationfile(string_noisefile)


if strcmp(string_noisefile, 'RGB-10-2-0.48-11111-64x32')
    params.pixelsize = 10;
    params.height = 32; params.width  = 64;
    params.refreshrate = 2;
    params.frames_pertrigger = 50;
    params.tstim = 2 * .0083275;
    params.type = 'RGB';
    params.RNG  = 11111;
end

if strcmp(string_noisefile, 'RGB-8-1-0.48-11111-80x40')
    params.pixelsize = 8;
    params.height = 40; params.width  = 80;
    params.refreshrate = 1;
    params.frames_pertrigger = 100;
    params.tstim             = .0083275;
    params.type = 'RGB';
    params.RNG  = 11111;
    
end


end