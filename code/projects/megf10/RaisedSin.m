function fit = RaisedSin(coef, X)

fit = coef(1) + coef(2)*(.5+.5*sind(X+coef(3))).^coef(4) ; % sin function

%Y = coef(3)+ coef(4) * exp(-((X-coef(1)).^2)/(2*coef(2))) ; % gaussian