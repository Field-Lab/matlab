function write_movie_idx(destination_mat,mov,mov_idx,stixel_sz)


dest=sprintf([destination_mat,'/%d'],mov_idx);
save([dest,'.mat'],'mov','-v7.3');
write_movie([dest,'.mat'],[dest,'.rawMovie'],stixel_sz);

end
