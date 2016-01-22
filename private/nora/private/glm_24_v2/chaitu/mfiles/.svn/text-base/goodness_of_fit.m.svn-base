% Goodness of fit of GLM parameters using time-rescaling and
% Kolmogorov-Smirnov Statistic

% If the model is good then the 'rescaled' ISI's obtained by by
% rescaling time by the computed CIF, should be i.i.d. exponentially
% distributed.

% Arguments:
% cifs for the data
% dt - time resolution of cifs
% sparse binary vector indicating spike times

function [uISIs uRISIs uLISIs KS nsp] = goodness_of_fit(cifs,kx,b,dt,D,plotFlag)

if (size(cifs,1) ~= size(D,1))
    fprintf('ERROR: cifs and spike train must be of the same size');
end

KS = zeros(size(cifs,2)); % Kolmogorov-Smirnov statistic on 1-D distribution test
nsp = zeros(size(cifs,2)); % number of samples for each neuron
1;
for i=1:size(cifs,2) %Iterate through neurons
    csum = cumsum(cifs(:,i)).*dt; % Cumulatively integrate the cif
    meanrate = sum(D(:,i))./(dt * size(D,1));
    
    kxmod = exp ( kx.*(log(meanrate)-b)/norm(kx) + b);
    
    kxsum = cumsum(kxmod).*dt; % Cumulatively integrate without the history terms - scale to match meanrate

    1;
    
    uISIs{i} = 1-exp(-sort(full(meanrate.*diff(find(D(:,i))).*dt),'ascend')); % Best poisson asssumption (rate = mean empirical rate)
    %uLISIs{i} = 1-exp(-sort(diff(kxsum(D(:,i))),'ascend')); % Best LNP assumption
    uRISIs{i} = 1-exp(-sort(diff(csum(D(:,i))),'ascend')); % Take the rescaled ISI's and sort them
    nsp(i) = length(uRISIs{i});
    
    yaxis = (0.5/nsp(i):1/nsp(i):((nsp(i)-0.5)/nsp(i)))';
    
    KS(i) = max(abs(uRISIs{i}-yaxis));
    if (plotFlag)
        %subplot(1,size(D,2),i);
        %plot([yaxis uISIs{i} uLISIs{i} uRISIs{i}],repmat(yaxis,1,4));
        plot([yaxis uISIs{i} uRISIs{i}],repmat(yaxis,1,3));
        %hold on, plot(uISIs{i},yaxis,'g--');
        %hold on, plot(yaxis,yaxis,'r-');
        title(sprintf('CDF of rescaled ISIs. Nsp=%d',nsp(i))), xlabel(sprintf('K-S test: %f/%f',KS(i),1.36/nsp(i)));
    end
    
end

uLISIs = {};
%figure, hist(diff(csum(D(:,i))));
%figure, plot(0:dt:(size(cifs,1)-1)*dt,csum)
    
    



