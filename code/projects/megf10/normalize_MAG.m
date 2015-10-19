function MAG_norm = normalize_MAG(DS)


MAG_norm = zeros(length(DS.num), length(DS.magmax));
MAG_all_temp = cell2mat(DS.MAG);
for cc = 1:length(DS.magmax)
    m = MAG_all_temp(:, cc)/max(MAG_all_temp(:, cc));
    MAG_norm(:, cc) = m;
end