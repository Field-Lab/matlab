
function filtered_lum2 = filter_lum(luminence_movie,filter_log)


filtered_lum2 = zeros(size(filter_log,2),size(luminence_movie,2));


for itau = 1: size(filter_log,2)
    filter_tau = filter_log(:,itau);
luminence_movie_ex = zeros(length(luminence_movie) + size(filter_log,1)-1,1);
luminence_movie_ex(length(filter_tau):end) = luminence_movie;

filtered_lum2(itau,:) = conv(luminence_movie_ex,filter_tau,'valid');
end

end