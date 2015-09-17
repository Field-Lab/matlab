function fval = f_val_population(Knew,Bnew,X,Y,dt,scale_cell) 
% calulate negative LL 
SU_inp_new = Knew'*X;
LL=0; T =size(X,2);
for icell = 1:size(Y,1)
    cell_inp = exp(SU_inp_new + repmat(Bnew(:,icell),[1,T]));
    cell_inp(:,sum(cell_inp,1)==0) = 0.0000001; % very small number to prevent it from  becoming NAN.
    LL = LL + (sum(cell_inp(:))*dt/T - (1/T)*log(sum(cell_inp,1))*Y(icell,:)')*T/scale_cell(icell);
end
fval = LL;
end