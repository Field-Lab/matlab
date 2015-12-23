function [ccf time]=correlation_shift(sp1, sp2, varargin) 
% shifted corrcoeff eg crosscorrelation function
%
%
%
% greschner

% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('dt', .004);
    p.addParamValue('offset', .15);
    p.addParamValue('shift', []);
    p.addParamValue('stop', []);
    p.addParamValue('start', []);
    p.addParamValue('max_isi', 0);
    p.parse(varargin{:});
    params = p.Results;
    
          
if ~isempty(params.start)
    sp1=sp1(find(sp1>=params.start));
    sp2=sp2(find(sp2>=params.start));
else
    params.start=min([sp1(1) sp2(1)]);
end    
if ~isempty(params.stop)
    sp1=sp1(find(sp1<=params.stop));
    sp2=sp2(find(sp2<=params.stop));
else
    params.stop=max([sp1(end) sp2(end)]);
end


if isempty(params.shift)
    params.shift=params.dt;
end
 

if params.max_isi~=0
    s=sp1;
    sp=zeros(size(s));
    sp(1)=1;
    sa=s(1);
    for i=2:length(s)
        if s(i)-sa>params.max_isi
            sp(i)=1;
            sa=s(i);
        end
    end
    sp1=s(find(sp));

    s=sp2;
    sp=zeros(size(s));
    sp(1)=1;
    sa=s(1);
    for i=2:length(s)
        if s(i)-sa>params.max_isi
            sp(i)=1;
            sa=s(i);
        end
    end
    sp2=s(find(sp));
end





if 0
    duration=max([max(spikes_1) max(spikes_2)]);

        spiketrain=zeros(duration/params.dt-params.offset/params.dt*2,params.offset/params.dt*2+1);
        t=histc(sp1,[0:params.dt:duration-params.dt]);
        spiketrain(:,1)=t(params.offset/params.dt+1:end-params.offset/params.dt);

        t=histc(sp2,[0:params.dt:duration-params.dt]);
        for i=1:params.offset/params.dt*2
            spiketrain(:,1+i)=t(i+1:end-params.offset/params.dt*2+i);    
        end
        temp=corrcoef(spiketrain);
        ccf=temp(1,2:end);

    time=[[-params.offset:params.dt:-params.dt] 0 [params.dt:params.dt:params.offset]];
end
if 0
    duration=datarun.duration/6;
    
        n=round(params.offset/params.dt);
        time=[[-params.dt*n:params.dt:-params.dt] 0 [params.dt:params.dt:params.dt*n]];
    
        tsp1=histc(sp1,[0:params.dt:duration-params.dt]);
        tsp2=histc(sp2,[0:params.dt:duration-params.dt]);
        
        for i=1:2*n
            t=[tsp1(n:end-n) tsp2(i:end-2*n+i)];
            temp=corrcoef(t);
            ccf(i)=temp(1,2);
        end
end

if 0 %binary histc  
    n=round(params.offset/params.dt);
    bins=[-n:n];
    time=bins*params.dt;


    a=histc(sp1,params.start:params.dt:params.stop-params.dt);
    t=find(a==1);
    if sum(t)/length(sp1)<99.999
        warning('correlation_shift: more than 0.1% of spikes were in bin with multiple entries');
    end
    a(a>1)=1;
    b=histc(sp2,params.start:params.dt:params.stop-params.dt);
    t=find(b==1);
    if sum(t)/length(sp2)<99.999
        warning('correlation_shift: more than 0.1% of spikes were in bin with multiple entries');
    end
    b(b>1)=1;
    
    ta=find(a);
    tb=find(b);
    lta=length(ta);
    ltb=length(tb);
    la=length(a);
    ma=lta/la;
    mb=ltb/la;
    
    %shift tb
    t=repmat(tb,1,length(bins))-repmat(bins,ltb,1);
    ltab=sum(ismember(t,ta)); 
    
    %cov(a,b)
    t=(1-ma)*(1-mb)*ltab + (1-ma)*-mb*(lta-ltab) + (1-mb)*-ma*(ltb-ltab) + -ma*-mb*(la-lta-ltb+ltab);
    
    %normalize->corrcoef
    ccf=t/(la-1)/sqrt(var(a)*var(b));    
end  


if 1 %fastest
  
    n=round(params.offset/params.shift);
    n=[-n:n];
    time=n*params.shift;
    
    bin=[params.start+time(end):params.dt:params.stop-time(end)];
    tsp1=histc(sp1,bin);
    
    for i=1:length(n)
        
        bin=[params.start+time(end)+time(i):params.dt:params.stop-time(end)+time(i)];
        tsp2=histc(sp2,bin);

        tt=corrcoef([tsp1 tsp2]);
        ccf(i)=tt(1,2);
    end  
   
end

if 0 % 50%slower
    tic
    n=round(params.offset/params.dt);
    n=[-n:n];
    time=n*params.shift;
    
    bin=[params.start:params.dt:params.stop];
    tsp1=histc(sp1,bin);
    tsp2=histc(sp2,bin);
    
    for i=1:length(n)
        tt=corrcoef([tsp1 circshift(tsp2,n(i))]);
        ccf(i)=tt(1,2);
    end  
    toc
end


if 0% much slower, not completed
    tic
    n=round(params.offset/params.shift);
    n=[-n:n];
    time=n*params.shift;
    
    bin=[params.start+time(end):params.dt:params.stop-time(end)];
    
    mat=zeros(length(bin),length(n)+1);
    mat(:,1)=histc(sp1,bin);
    mat(:,2:end)=repmat(histc(sp2,bin),1,length(n));
    
    mat=circshift(mat,[0 n]);
    
    t=corrcoef(mat);
    
    ccf=t(2:end,1);

    toc
end
    





    

















    
    