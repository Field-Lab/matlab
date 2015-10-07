
function[ predictedResponse,lam] = predictGMLM_full_gains(fitGMLM,mov_use,nTrials,nl,g2)
filters = fitGMLM.Linear.filter;
nFrontEnds = length(filters);
mu=fitGMLM.mu;
hbas= fitGMLM.hist.hbasis;
h_filt = fitGMLM.hist.hBas;

mov_use_big = zeros(size(mov_use,1),10*size(mov_use,2));

for istimdim=1:size(mov_use,1)
xx=mov_use(istimdim,:);
xx=repmat(xx,[10,1]);
xx=xx(:)';
mov_use_big(istimdim,:)=xx;
end


if(length(filters{1}) > size(mov_use_big,1)) % Add extra 1 in end.
    mov_use_big = [mov_use_big',ones(size(mov_use_big,2),1)]';
end

mov_filtered=mov_use_big;
kx=cell(nFrontEnds,1);

for ifilter=1:nFrontEnds
    kx{ifilter} = exp(filters{ifilter}'*mov_filtered);
end



lam=kx{1};
for ik=2:length(kx)
lam=lam+kx{ik};
end
lam=lam+mu;

lam=nl(g2*lam)*(1/1200);
predictedResponse = zeros(nTrials,length(lam));

% for itrial=1:nTrials
%     itrial
%     for itime=1:size(predictedResponse,2)
%         
%         yhistory=zeros(size(hbas,1),1);
%         if(itime>size(predictedResponse,2))
%         yhistory = predictedResponse(itrial,itime-1:time-size(hbas,1))';
%         end
% predictedResponse(itrial,itime)=poissrnd(lam(itime)) * exp(h_filt'*hbas'*yhistory);    
%     end
% end
% 
% predictedResponse= predictedResponse';

nbins=length(lam);
PS = fitGMLM.hist.hexpanded;
  cif_psgain = exp(PS);
    ps_bins     = length(cif_psgain);
    for i_trial = 1 : nTrials
       
        cif_ps       = lam;
        binary_simulation = zeros(1,nbins);
        for i = 1 : nbins - ps_bins;
            roll = rand(1);
            if roll <(cif_ps(i));
                cif_ps(i+1: i + ps_bins) =  cif_ps(i+1: i + ps_bins) .* (cif_psgain');
                binary_simulation(i)= 1;
            end
        end
       predictedResponse(i_trial,:) = binary_simulation ;
    end
 
predictedResponse= predictedResponse';
end