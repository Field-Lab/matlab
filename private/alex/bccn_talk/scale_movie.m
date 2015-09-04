function scaled_movie = scale_movie(raw_movie, stix_size)

scaled_movie=zeros((320/stix_size)^2, 3600);

stix_size = stix_size/2;
% average in space

cnt=1;
for i=1:stix_size:160
    for j=1:stix_size:160
        tmp=reshape(raw_movie(i:i+stix_size-1,j:j+stix_size-1, :), stix_size*stix_size, 3600);
        scaled_movie(cnt,:)=mean(tmp);
        cnt=cnt+1;
    end
end

% % average in time
% scaled_movie = zeros((320/stix_size)^2, 3600/frames_to_average);
% cnt = 1;
% for i=1:frames_to_average:3600
%     scaled_movie(:, cnt) = mean(movie_tmp(:,i:i+frames_to_average-1),2);
%     cnt = cnt+1;
% end
