OutputPath = 'Z:\Szarowaty\Projekt_neuro';
data_id = 'PD_data_2014_05_15_16_10';

events = reshape(load([OutputPath filesep data_id filesep data_id '.dat']),[],2);

figure(1)
x = -945:10:945;
y = -450:10:450;
for e = 1:length(events)
    load([OutputPath filesep data_id filesep 'PD_p' num2str(events(e,2)) '_m' num2str(events(e,1)) '.mat']);
    if ~orient
        y = polyval(line,x);
    else
        x = polyval(line,y);
    end
    plot(x,y);
    title(num2str(events(e,2)))
    hold on
    axis([-945 945 -450 450])
    %pause
end