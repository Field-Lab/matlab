DataPath = 'Z:/Data/scan_new';
PatternNumber = [1:64,321:384]; 
BadElectrodes = [353 354 420];
%MovieNumber = 51:2:63;
MovieNumber = 63;
Resolution = 10;

sigma = zeros(length(MovieNumber),length(PatternNumber),1); % last dimension will be modified in loop
amplitude = zeros(length(MovieNumber),length(PatternNumber),1);
offset = zeros(length(MovieNumber),length(PatternNumber),1);
r_mean = zeros(length(MovieNumber),length(PatternNumber),1);
dists = zeros(length(MovieNumber),length(PatternNumber),1);

fx = @(a,x_0,s,y_0,x)(y_0 + a*exp(-(x-x_0).^2./(2*s.^2)));

figure(1)
for m = 1:length(MovieNumber)
   for p = 1:length(PatternNumber)
       n = 1;
       [plots,Distances] = getCompleteIntersectionSet (DataPath,PatternNumber(p),MovieNumber(m),Resolution,BadElectrodes);
       for d = 1:length(Distances)
           pp = squeeze(plots(d,:,:));
           pp(isnan(pp(:,2)),:) = [];
           if length(pp(:,2)) > 1
               f = fit(pp(:,1),pp(:,2),fx,'StartPoint',[400,0,200,0.1]);
               sigma(m,p,n) = f.s;
               amplitude(m,p,n) = f.a;
               offset(m,p,n) = f.y_0;
               r_mean(m,p,n) = f.x_0;
               dists(m,p,n) = Distances(d);
               n = n + 1;
           end
       end
   end
   plot(reshape(dists(m,:,:),[],1),reshape(sigma(m,:,:),[],1),'.')
   xlabel('perpendicular distance [um]')
   ylabel('sigma [um]')
   title(strcat('Sigma for MovieNumber = ',num2str(MovieNumber(m))));
   file_sigma = strcat('output/m',num2str(MovieNumber(m)),'_sigma.png');
   print('-dpng',file_sigma,'-r120');
end