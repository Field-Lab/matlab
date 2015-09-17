function [fval,grad] = minL_dec3D(x,K,B,Y,dt,d,T,scale_cells,ttf,Tker)

x =reshape(x,[d,T]);
 x= filterMov_cone(x,logical(ones(d,1)),squeeze(ttf));
      
%% objective value

SU_inp = K'*x;
fval = 0;
for icell = 1:size(Y,1)
    cell_su_inp = exp(SU_inp + repmat(B(:,icell),[1,T]));
fval = fval+(sum(cell_su_inp(:))*dt/T - Y(icell,:)*log(sum(cell_su_inp,1))'/T)*(T/scale_cells(icell));
end
fval = fval;

%% grad
% 
Ns=size(K,2);

grad = zeros(size(K,1),T);
for icell=1:size(Y,1)
    cell_su_inp = exp(SU_inp + repmat(B(:,icell),[1,T]));
    grad = grad + (K*cell_su_inp*dt/T - K*((cell_su_inp./repmat(sum(cell_su_inp,1),[Ns,1])).*repmat(Y(icell,:),[Ns,1]))/T)*T/scale_cells(icell);
end


grad = grad*Tker;

grad = grad(:) ;
% grad=zeros(d,T);
% for icell=1:size(Y,1)
%     cell_su_inp = exp(SU_inp + repmat(B(:,icell),[1,T]));
%     
% end
% grad=grad(:);

grad=gather(grad);
fval= gather(fval);
end