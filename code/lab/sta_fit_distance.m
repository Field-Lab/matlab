function [dist normdist] = sta_fit_distance(sta_fit_a,sta_fit_b)
%distance and normalized distance of fits
%
%greschner

if isempty(sta_fit_a)|isempty(sta_fit_b)
    warning('sta_fit empty')
    dist=NaN; 
    normdist=NaN;
else
    m=sta_fit_a.mean-sta_fit_b.mean;
    
    [t,rho]=cart2pol(m(1),m(2));
    sig_a=sta_fit_a.sd(1)*sta_fit_a.sd(2)./sqrt(sta_fit_a.sd(1)^2*sin(t+sta_fit_a.angle).^2+sta_fit_a.sd(2)^2*cos(t+sta_fit_a.angle).^2);

    [tt,rho]=cart2pol(-m(1),-m(2));
    sig_b=sta_fit_b.sd(1)*sta_fit_b.sd(2)./sqrt(sta_fit_b.sd(1)^2*sin(t+sta_fit_b.angle).^2+sta_fit_b.sd(2)^2*cos(t+sta_fit_b.angle).^2);

    dist=rho;
    
    normdist=2*rho/(sig_a+sig_b);
end
 
    
    

%test plot
if 0
    clf
    hold on

    plot_rf_outlines_(sta_fit_a);
    plot_rf_outlines_(sta_fit_b);
    
    set(gca,'DataAspectRatio',[1 1 1]);

    plot(sta_fit_a.mean(1),sta_fit_a.mean(2),'.r');
    plot(sta_fit_b.mean(1),sta_fit_b.mean(2),'.r');
    
    [xa,ya]=pol2cart(t,sig_a);
    [xb,yb]=pol2cart(tt,sig_b);
    
    plot([xa+sta_fit_a.mean(1) xb+sta_fit_b.mean(1)],[ya+sta_fit_a.mean(2) yb+sta_fit_b.mean(2)]);
    
    plot(xa+sta_fit_a.mean(1),ya+sta_fit_a.mean(2),'.g');
    plot(xb+sta_fit_b.mean(1),yb+sta_fit_b.mean(2),'.g');

    [xa,ya]=pol2cart(t+pi,sig_a);
    [xb,yb]=pol2cart(tt+pi,sig_b);
    
    plot([xa+sta_fit_a.mean(1) xb+sta_fit_b.mean(1)],[ya+sta_fit_a.mean(2) yb+sta_fit_b.mean(2)]);
    
    plot(xa+sta_fit_a.mean(1),ya+sta_fit_a.mean(2),'.k');
    plot(xb+sta_fit_b.mean(1),yb+sta_fit_b.mean(2),'.k');
    
end

 
   