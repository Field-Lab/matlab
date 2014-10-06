%%
for i = 1:20
    Cs{1,i} = squareconnect(i);
    Cs{2,i} = hexconnect(i);
    Cs{3,i} = squarediagconnect(i);
end; clear i
cs = [4 6 8];

%%
Bs = [0.01 0.1 0.2 0.5 1 1.5 2 2.5 3 5 10 100 1000];
save netopt

%%
for ci = 1:3
    c = cs(ci);
    for l = 1:20
        for bi = 1:length(Bs)
            b = Bs(bi);
            disp([num2str(c) ',' num2str(l) ',' num2str(b)]);
            
            W = cellnet(Cs{ci,l}, b);
            Ns(ci,l,bi) = 1 / sum(W(:,1).^2);
            Nsovercs(ci,l,bi) = N(ci,l,bi) / c;
        end; clear bi b W
    end; clear l
end; clear ci c
save netopt