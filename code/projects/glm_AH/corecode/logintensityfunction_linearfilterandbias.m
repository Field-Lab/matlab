function [lcif_mu , lcif_kx ] = logintensityfunction_linearfilterandbias (pstar , Paramind , Movie, spikebins_perframe, LinearFilter) 
%CHECKED 2012-12-06    AK HEITMAN

%INPUTS:
% pstar    = list of final parameter values
% Paramind.  a hack to parse pstar into the correct ones
%            needs to have .L , .MU for this file
% Movie   movie file (xyspace, time)    
% LinearFilter (xyspace,time )   (only if it has already been computed)
% spikebins_perframe : no spiking, but want the lcifs to be of correct
% dimension for use in simulatoin     >=1 and integer valued!

% Truly just a wrapper for filterstimulus_train3AH



% mu is just a constant and needs to be put into vector for each spike bin
mu     = pstar ( Paramind.MU);
frames = size(Movie,2);
mu_bin = mu * ones (frames * spikebins_perframe,1); 

% kx is the stimulus run through linear filter , at stimframe rate
if exist ('LinearFilter' , 'var')
    string_mode = 'fixed_filter';
	kx          = filterstimulus_train3AH(pstar,Paramind,string_mode,Movie,LinearFilter);
end
if ~exist ('LinearFilter' , 'var')
    string_mode = 'rk2';
	kx          = filterstimulus_train3AH(pstar,Paramind,string_mode,Movie);
end

% put into spike bins sized time
kx_bin = zeros(frames*spikebins_perframe , 1);
sbpf   = spikebins_perframe;
for i = 1 : length(kx)
        indices = 1 + (i-1)*sbpf: i*sbpf ;
        kx_bin(indices) = kx(i);
end

% final assignment of output
lcif_mu = mu_bin;
lcif_kx = kx_bin;

end




    



