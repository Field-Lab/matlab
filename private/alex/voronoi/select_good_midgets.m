datarun = load_data('/Volumes/Analysis/2010-09-24-1/data006/data006');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

cone_template = 0;
border = 10;
for visionID = datarun.cell_types{4}.cell_ids
    datarunID = find(datarun.cell_ids == visionID);
    sta = squeeze(datarun.stas.stas{datarunID}(:,:,1,:));
    sta = sta(:,:,4);
    [a,b] = find(sta == min(sta(:)));
    if a>border && a<320-border && b>border && b<320-border
        cone_template = cone_template + sta(a-border:a+border, b-border:b+border);
    end
end

figure
colormap gray
imagesc(cone_template)

figure
plot(cone_template(11,:))
hold on
plot(cone_template(:,11))