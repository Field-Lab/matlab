function population = map_population_electrodes(population,elecs)

nElecs = length(elecs.x);
population.elecs.dist_elecs = zeros(population.nCones,nElecs);
population.elecs.weight_elecs = zeros(population.nCones,nElecs);
for icell=1:population.nCones
    xCell =population.conesX;
    yCell =population.conesY;
    
    dist_elecs = (elecs.x - xCell(icell)).^2+(elecs.y - yCell(icell)).^2;
    population.elecs.dist_elecs(icell,:) = dist_elecs;
    population.elecs.weight_elecs(icell,:) =   population.elecs.dist_elecs(icell,:)<80;
    
end


end