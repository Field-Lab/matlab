function tb=trapez_blackman(length,marg);

m=blackman_ph(marg*2)';
size(m)
tb=ones(1,length);
tb(1,1:marg)=m(1,1:marg);
tb(1,(length-marg+1):length)=m(1,(marg+1):marg*2);

