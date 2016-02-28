function population = map_population_electrodes_subunits(population,elecs)

%% give weights to elec SU
nSU = elecs.nSU;
population.elecs.dist_elec_su = zeros(population.nCones,nSU);
population.elecs.weight_elec_su = zeros(population.nCones,nSU);
for icell=1:population.nCones
    xCell =population.conesX;
    yCell =population.conesY;
    
    dist_elec_su = (elecs.suX - xCell(icell)).^2+(elecs.suY - yCell(icell)).^2;
    population.elecs.dist_elec_su(icell,:) = dist_elec_su';
    [v,iidx] = sort(dist_elec_su,'ascend');
    xx = zeros(1,nSU);
    % xx(iidx(1:6))=1;
    xx(dist_elec_su<20)=1;
    population.elecs.weight_elec_su(icell,:) =   xx;
    
end

%% figure? 
figure;
plot(population.conesX,population.conesY,'r.');
hold on;
plot(elecs.x,elecs.y,'b.');
hold on;
plot(elecs.suX,elecs.suY,'g.');
hold on;
for isu=1:elecs.nSU
iidx = elecs.su_elec(isu,:);
xx = elecs.x(logical(iidx));
yy = elecs.y(logical(iidx));
aa = convhull([xx,yy]);
plot(xx(aa),yy(aa),'Color',[0.8,0.8,0.8]);
hold on;
end

iidx=1:elecs.nSU;
cols = distinguishable_colors(population.nCones);
for icell=1:population.nCones
hold on;
sus = iidx(logical(population.elecs.weight_elec_su(icell,:)));
for isu=1:length(sus)
plot([population.conesX(icell),elecs.suX(sus(isu))],[population.conesY(icell),elecs.suY(sus(isu))],'Color',cols(icell,:));
hold on;
end

end

end