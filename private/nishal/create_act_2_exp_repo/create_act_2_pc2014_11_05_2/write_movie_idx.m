function write_movie_idx(destination_mat,mov,mov_idx)


dest=sprintf([destination_mat,'/%d'],mov_idx);
save([dest,'.mat'],'mov');
write_movie([dest,'.mat'],[dest,'.rawMovie'],10);

end
