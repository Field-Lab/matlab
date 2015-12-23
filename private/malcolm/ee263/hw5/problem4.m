% problem 4

x0 = [Tw; Ta*ones(n,1)];
z0 = x0 - Ta*ones(n+1,1);

% given P, calculate temperature of espresso at P+15 seconds
P = 5:0.1:15;
Tfinal = nan(length(P),1);
for i = 1:length(P)
    zP = expm(P(i)*A)*z0; % wait P seconds
    zP(1) = Te-Ta; % water exchanged for espresso
    zP15 = expm(15*A)*zP; % wait 15 seconds
    Tfinal(i) = zP15(1)+Ta;
end

Popt = P(Tfinal==max(Tfinal));
Topt = max(Tfinal);