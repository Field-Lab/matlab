function plotio_curve_meancorrected(aa,resp,pts,colr)
bin=[];

 for iprc=100/pts:100/pts:100
 bin=[bin; prctile(aa,min(iprc,100))];
 end
 
 meanR=[];meanaa=[];varaa=[];
 for ibin=1:length(bin)
 if(ibin==1)
     tms = aa<bin(ibin);
 else
     tms= aa>bin(ibin-1) & aa<bin(ibin);
 end

 meanR = [meanR ;mean(resp(tms))];
 meanaa  = [meanaa;mean(aa(tms))];
varaa=[varaa;sqrt(var(aa(tms))/numel(aa(tms)))];
 end
 meanR=meanR - mean(resp);

 errorbar(meanaa,meanR,varaa,colr)
%plot(meanaa,meanR,colr)
 
end