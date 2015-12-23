function w=prox_binary(v,lam,u,l)
% u=127.5
% l=-127.5
% lam =1.2?

w=0*v;
w(v>=u+2*lam)=v(v>=u+2*lam)-2*lam;
w(v<=l-2*lam)=v(v<=l-2*lam)+2*lam;
w(v<=u/lam & v>=l/lam)=v(v<=u/lam & v>=l/lam)*lam;
w(v>u/lam & v<u+2*lam) = u;
w(v<l/lam & v>l-2*lam) = l;

end