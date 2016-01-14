function save_parfor_stas(filename, sta, e_sta)
%SAVE_PARFOR_STAS Saves stas inside a parfor loop

save(filename, 'sta', 'e_sta');

end % save_parfor_stas