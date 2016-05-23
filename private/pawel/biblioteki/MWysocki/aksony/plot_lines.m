file{1} = 'Z:\Szarowaty\Projekt_neuro\PD_retina1_2014_05_16_14_33_m63.mat';
file{2} = 'Z:\Szarowaty\Projekt_neuro\PD_retina1_2014_05_16_10_42_m51.mat';
file{3} = 'Z:\Szarowaty\Projekt_neuro\PD_retina1_2014_05_16_14_52_m57.mat';

x = -945:10:945;
y = -450:10:450;

figure(1)
hold off
for f = 1:length(file)
    load(file{f})
    for p = 1:length(patterns)
        x = -945:10:945;
        y = -450:10:450;
        if ~isempty(rms{p}) && rms{p} < 10
            if orient{p}
                y = polyval(line{p},x);
            else
                x = polyval(line{p},y);
            end
            plot(x,y);
            hold on
            axis([-945 945 -450 450])
        end
    end
end