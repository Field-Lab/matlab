%% avg RF
function avg_profile_NB(datarun, cell_spec)

cells = get_cell_indices(datarun, cell_spec);

dataruns{1} = load_sta(dataruns{1});
avg_rf = sum(get_average_rf(dataruns{1}, cell_type, 'scale', 5),3);

n_angles = 50;
slice = zeros(200,1);
for i = 1:n_angles
    angle = i*360/n_angles;
    temp = imrotate(avg_rf, angle, 'crop');
    slice = slice+ improfile(temp, [200 200], [1 200]);
end
slice = slice/n_angles;
a = fit((1:200)'/5,slice, fittype('gauss1'));
width = a.c1/sqrt(2);
avg_profile = [slice, (1:200)'/(5*width)];

end