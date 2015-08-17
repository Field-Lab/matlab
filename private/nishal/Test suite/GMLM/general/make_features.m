function [features,filter_log,filtered_lum] = make_features(maskedMovdd_sliced,movie_full_sliced)

[filtered_lum,filter_log,times,luminence_movie] = luminence_signal2(movie_full_sliced+0.5+0.01);
 
 
 features = maskedMovdd_sliced;
 xdim = size(maskedMovdd_sliced,1);
 features = [features;maskedMovdd_sliced./repmat(filtered_lum,[xdim,1])];
 %features = [features;maskedMovdd_sliced./repmat(filtered_lum.^2,[xdim,1])];

 features = [features;maskedMovdd_sliced.*repmat(filtered_lum,[xdim,1])];
 %features = [features;maskedMovdd_sliced.*repmat(filtered_lum.^2,[xdim,1])];

 
end