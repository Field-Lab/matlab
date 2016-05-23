dataPath = 'D:\Home\Pawel\analysis\retina\2012-09-27-4\analysis_2014_05_16\MWysocki';
dataName = 'retina1_2014_05_16_17_31'; % ZMIENIC na odpowiednia date
                    % (dataName jest po prostu nazwa dowolnego z plikow z danymi,
                    % ale po odcieciu 'PD_' z przodu i '_m<numer_moviesa>.mat' z tylu)
movieRange = 51:2:63;   % chyba dla takich moviesow zbieralismy dane

x = -945:10:945;
y = -450:10:450;

figure(5)
hold off
for iMovie = 1%1:length(movieRange)
    file = [dataPath filesep 'PD_' dataName '_m' num2str(movieRange(iMovie)) '.mat']
    load(file);
    % Do workspace'a trafiaja teraz zmienne: 'line','orient','v0','ve','vp','edges','func','success','rms','patterns'
    % patterns to tablica, a reszta jest klasy CellArray
    for p = 1:length(patterns)
        x = -945:10:945;
        y = -450:10:450;
        if ~isempty(rms{p}) && rms{p} < 10 % kryterium odrzucania nieudanych lub 'brzydkich' przypadkow (zamiast 10 um mozna wstawic inna wartosc)
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