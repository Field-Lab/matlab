function info = compute_autocor_renewal(spikearray,TP,K)
%*************
%****** info = compute_autocor_renewal(spikearray,timediv,K)
%*** inputs:
%****** spikearray - spike matrix (trials x timepoints)
%****** TP - the discretization of time, 0.5 ms
%****** K - number of basis to use smoothing kernel
%******
%****** NOTE:  compute Poisson shuffled trial-trial versions
%******        that preserve the rate across time in trials
%******        and subtract out that part of correlation
%******        before packaging (correct for non-stationary rate)
%******
%***** outputs:
%***
%    info.timediv = TP * (1:M2);   % Time delays of Autocorrelation
%    info.autocor = Autocorrelation, Prob(spike(t)|spike(0))
%    info.autocor_s = Jacknife estimate of confidence intervals
%    info.jackcor = examples from Jacknifing across trials
%    info.jackrate = mean rate for given Jacknife 
%    info.poisscor = examples from Poisson tests
%    info.renewcor = examples from Renewal tests
%    info.renewrate = rate of each example 
%    info.pautocor = Mean Autocorrelation of several rate matched Poisson
%    info.sautocor = Std of several rate matched Poissons
%    info.rate = height of Poisson = mean rate
%*****
%********** added features for renewal process
%    info.rautocor = Mean Autocorrelation of several renewal processes
%                     (spike trains of same ISI's but order shuffled)
%    info.srautocor = Std of several matched renewal processes
%************************************

  info = []; 
  N = size(spikearray,1);  % number of trials
  NT = 40;  % estimate poisson rate matched distributions
             % from random shuffles of original spike trains
            %*** and also estimated ISI shuffle renewal autocor
  M = floor( size(spikearray,2)*3/4);  %make maximum ISI delay 
                                      % 3/4 of total interval
  M2 = M;
  pstart = 1;
  pend = size(spikearray,2);
   
  %****************************** for norm scale, use 0.5ms divisions
   
  %**************** find the autocor distributions ********************
  %********* (also, correct for there being less samples of longer
  %*********  ISI durations due to finite size of data interval by
  %*********  dividing by the total number of possible samples at
  %*********  each ISI, for example, for 800ms, you can sample a
  %*********  1ms ISI 800 times, but a 400ms isi only 400 times).
  iisicnt = zeros(1,M2);
  iisibin = zeros(1,M2);
  %*********
  zisicnt = zeros(1+(2*NT),M2);
  zisibin = zeros(1+(2*NT),M2);
  szisibin = zeros(1+(2*NT),M2);
  %************ get Jacknife estimate of confidence intervals actual data
  zisicnt_trial = zeros(N,M2);  % estimate per trial to get Jacknife
  zisibin_trial = zeros(N,M2);
  zisijack = zeros(N,M2);  % jacknife estimates (leaving one trial out at a time)
  szisijack = zeros(N,M2); % smoothed isi jacknife
  renewrate = [];
  
  for nkk = 1:(1+(2*NT))
        
     if (nkk>1)  
         if (nkk <= (1+NT))  % rate matched poisson
           % estimated one of NT samples of Poisson rate matched spike
           %******* array by random shuffling of existing spike array
           anospikes = zeros(size(spikearray));
           
           % shift-predictor ........
           %********* first, jitter times at least 10 ms as
           %********* you only want to remove low freq trends
             MZ = size(spikearray,2);
             sparray = zeros(size(spikearray));
             if (1)
               sparray = spikearray;
             else
               for kz = 1:size(spikearray,1)
                  y = find( spikearray(kz,:) > 0);
                  for kt = 1:size(y,2)
                      aa = max(1,(y(kt)-10));
                      bb = min(MZ,(y(kt)+10));
                      rep = 1;
                      while (rep == 1)
                         it = aa + floor( (bb-aa) * rand);
                         if (sparray(kz,it) == 0)
                             sparray(kz,it) = 1;
                             rep = 0;
                         end
                      end
                  end
               end
             end
             %********* compute low freq due to rate variation      
             for ii = 1:size(spikearray,2)
                rr = randperm(size(spikearray,1));
                anospikes(:,ii) = sparray(rr,ii);
             end
             %**********************************************
             
           
         else  % compute autocorrelation of renewal process
             
           anospikes = shuffleisi(spikearray);  % random renewal process
           renewrate = [renewrate mean(mean(anospikes))];  % mean rate
         end
         
     else  % and otherwise, use actual array itself
          anospikes = spikearray;   %use original spike array first time
     end
     
     %**********************************************
     for ii = 1:N  %ii is the trial number
       
       %********** report progress during computation ********  
       if (mod(ii,20) == 0)
         if (nkk == 1)
             disp(sprintf('Computing autocor, trial: %d',ii));
         else
           if (nkk <= (1+NT))  
             disp(sprintf('Computing Rate Matched Poisson %d, test: %d',(nkk-1),ii));    
           else
             disp(sprintf('Computing Renewal of Process %d, test:%d',(nkk-1-NT),ii));
           end
         end
       end
       %***************************************************
       
       jj = pstart;
       while (jj<pend)
           if (anospikes(ii,jj)==1)
               firsthit = 0;
               for kk = (jj+1):pend
                 if ((kk-jj)<=M)
                    ikk = (kk-jj);
                    %********************************************* 
                    if (anospikes(ii,kk)==1)
                        zisibin(nkk,ikk) = zisibin(nkk,ikk) + 1;
                        if (nkk == 1)
                           zisibin_trial(ii,ikk) = zisibin_trial(ii,ikk) + 1;
                           if (firsthit == 0)
                               firsthit = 1;
                               iisibin(ikk) = iisibin(ikk) + 1;
                           end
                        end
                    end
                    zisicnt(nkk,ikk) = zisicnt(nkk,ikk) + 1;
                    if (nkk == 1)
                        zisicnt_trial(ii,ikk) = zisicnt_trial(ii,ikk) + 1;
                        iisicnt(ikk) = iisicnt(ikk) + 1;
                    end
                 end
               end              
           end
           jj = jj + 1;
       end %jj updated for intervals over a single trial
     
     end  % loop over all trials in matrix
     
     %************ normalize each bin by # of samples at given ISI
     %****** gives probability of spike at ISI given spike at time 0
     for ii = 1:size(zisibin,2)
         if (zisicnt(nkk,ii))
            zisibin(nkk,ii) = zisibin(nkk,ii) / zisicnt(nkk,ii);
         else
            zisibin(nkk,ii) = 0;
         end
     end
         
     %************ now for nkk == 1 (original data) build up the
     %************ Jacknife estimate of confidence bounds ******
     if (nkk == 1)

        for ii = 1:N
         xx = find( (1:N) ~= ii );
         sumo = sum( zisibin_trial(xx,:) );
         cumo = sum( zisicnt_trial(xx,:) );
         zisijack(ii,:) = sumo ./ cumo;
         zjackrate(ii) = mean(mean( spikearray(xx,:) ));
        end
     end
     
 end % loop of nkk, 1 original train, others Poisson random samples
    
 
 %********* compute Poisson mean estimate ***************
 info.poissest = mean(zisibin(2:(NT+1),:));      
 info.rate = mean(mean( spikearray ));
 info.pdeviate = info.poissest - info.rate;  % deviation from Poisson due to within
                                             % trial drift in the firing 
 %******** now subtract out that component from everything ******
 %******** and if desired, do smoothing on raw auto once corrected
 for nkk = 1:(1+(2*NT))   % for autocor and poisson
                      % not for renewal process, which by def, has no 
                      % trend in mean rate over trial
     if (nkk <= (1+NT) )                 
       zisibin(nkk,:) = zisibin(nkk,:) - info.pdeviate;
     end
     %*********** compute a smoothed version *************
     if (nkk == 1)
         [powo,basis] = compute_logtime_smooth(zisibin(nkk,:),K,[]);
         szisibin(nkk,:) = powo;
     else
         szisibin(nkk,:) = compute_logtime_smooth(zisibin(nkk,:),K,basis);
     end
 end
 for i = 1:N
     zisijack(i,:) = zisijack(i,:) - info.pdeviate;
     szisijack(i,:) = compute_logtime_smooth(zisijack(i,:),K,basis);
 end
 autocor_s = sqrt( (N-1) * var( zisijack ) );   % Jacknife estimate
 zautocor_s = sqrt( (N-1) * var( szisijack ) );  % Jacknife of smoothed
 
 %**********************************************
 info.timediv = TP * (1:M2);   % Time delays of Autocorrelation
 info.autocor = zisibin(1,:);  % Autocorrelation, Prob(spike(t)|spike(0))
 info.autocor_s = autocor_s;
 %********
 info.poisscor = zisibin(2:(NT+1),:);
 info.jackcor = zisijack;
 info.jackrate = zjackrate;
 info.renewcor = zisibin((NT+2):(1+(2*NT)),:);
 info.renewrate = renewrate;
 %********
 info.pautocor = mean( zisibin(2:(NT+1),:) );  % mean Poisson
 info.sautocor = std( zisibin(2:(NT+1),:) );   % std of Poissons
 info.rautocor = mean( zisibin((NT+2):(1+(2*NT)),:) );  % mean renewal
 info.srautocor = std( zisibin((NT+2):(1+(2*NT)),:) );   % std of renewal
 
 %******************** smoothed versions of same thing 
 info.zautocor = szisibin(1,:);  % Smooth Autocorrelation, Prob(spike(t)|spike(0))
 info.zautocor_s = zautocor_s;
 info.zpautocor = mean( szisibin(2:(NT+1),:) );  % mean Poisson
 info.zsautocor = std( szisibin(2:(NT+1),:) );   % std of Poissons
 info.zrautocor = mean( szisibin((NT+2):(1+(2*NT)),:) );  % mean renewal
 info.zsrautocor = std( szisibin((NT+2):(1+(2*NT)),:) );   % std of renewal
  
 %***************** send back ISI distribution for reference
 info.isidist = iisibin ./ iisicnt;
 info.isidist = info.isidist / sum(info.isidist);
 
return;

 
%********************************************
function anospikes = shuffleisi(allspikes)
%******* function anospikes = shuffleisi(allspikes)
%**** creates a renewal process spike train from original train
%*** input: allspikes (trials x time)  0s and 1s
%*** output: identical sized spike array, ISI shuffled each trial

  anospikes = zeros(size(allspikes));
  N = size(allspikes,1);
  TT = size(allspikes,2);
  
  for ni = 1:N % trials
      
           %******* first, collect all interspike intervals
           inters = [];
           lastspike = 0;
           for tt = 1:TT
               if (allspikes(ni,tt)>0)
                   inters = [inters (tt-lastspike)];
                   lastspike = tt;
               end
           end
           if (lastspike ~= TT)
               inters = [inters (TT-lastspike)];
           end
           
           %******** then make a random permutation ******
           mo = size(inters,2);
           rr = randperm(mo);
           
           %********** and write back in different order
           lastspike = 0;
           for tt = 1:(mo-1)
               lastspike = inters(rr(tt)) + lastspike;
               anospikes(ni,lastspike) = 1;
           end
  end
  
return;
