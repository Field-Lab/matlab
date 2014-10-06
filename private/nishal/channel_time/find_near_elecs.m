function nearest_elecs = find_near_elecs(positions,elec)
relative_pos=positions-repmat(positions(elec,:),size(positions,1),1);

dist=sum(relative_pos.^2,2);

[sort_dist,I]=sort(dist);

nearest_elecs=I(1:7);

end