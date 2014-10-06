function [Nc] = convert_matrix(N)    

    NcRows = 0;    
    for (i = 1:size(N, 1))
        NcRows = NcRows + N(i, 4);   
    end
    Nc = zeros(NcRows, 3);    %tworze macierz o odpowiedniej wielkosci=l.czasow
    
    timeC = 1;
    for (i = 1:size(N, 1))
        if (i ~= 1)
            if (N(i, 1) ~= N(i - 1, 1))
                timeC = timeC + 1;        
            end
        end
    end                       %przejscie do kolejnych, roznych od poprzednich czasow
    
    timeA = zeros(timeC, 1);
    timeC = 1;
    for (i = 1:size(N, 1))
        if (i ~= 1)
            if (N(i, 1) ~= N(i - 1, 1)) 
                timeA(timeC) = i - 1;   %wiersze, gdzie konczy sie dany czas
                timeC = timeC + 1;
            end
        end
    end
    timeA(timeC) = i;          
    
    NcIndex = 1;
    for (i = 1:timeC)
        for (j = 1:timeA(i))       %przez wszystkie wiersze, po kolei od 1 do 4, od 1 do(tam gdzie sie konczy kolejny czas)
            if (N(j, 4) > 0)
                Nc(NcIndex, 1) = N(timeA(i), 1);     
                Nc(NcIndex, 2) = N(j, 2);
                Nc(NcIndex, 3) = N(j, 3);                 %powielam kolejne wiersze tyle razy, ile maja w 4 kolumnie
                N(j, 4) = N(j, 4) - 1;
                NcIndex = NcIndex + 1;
            end
        end
    end    
end