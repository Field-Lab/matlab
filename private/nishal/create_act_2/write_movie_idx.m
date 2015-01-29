function write_movie_idx(destination_mat,mov,mov_idx,stixel_sz)


dest=sprintf([destination_mat,'/%d'],mov_idx);
save([dest,'.mat'],'mov');
write_movie([dest,'.mat'],[dest,'.rawMovie'],stixel_sz);

end
