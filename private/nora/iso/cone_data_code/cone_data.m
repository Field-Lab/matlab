rod = table2array(readtable('/Users/Nora/Desktop/mf-lmr_1', 'Delimiter', '\t', 'ReadVariableNames', false));
cone = table2array(readtable('/Users/Nora/Desktop/mf-lms_1', 'Delimiter', '\t', 'ReadVariableNames', false));

L = cone(:,1);
M = cone(:,2);
S = cone(:,3);
R = rod(:,3);

photoreceptor.spectra = [L M S R];
photoreceptor.column_labels = 'columns are L M S Rod';
photoreceptor.row_labels = 'rows are 370:730 nm';
photoreceptor.source = 'Baylor, Nunn and Schnapf 1987';
photoreceptor.species = 'macaque';
save('macaque_photoreceptor_spectra.mat', 'photoreceptor')