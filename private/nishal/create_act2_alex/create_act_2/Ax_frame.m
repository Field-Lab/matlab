function resp=Ax_frame(stas,frame)

itimeSta=1;

n_cell=length(stas);
resp=zeros(n_cell,1);

parfor icell=1:n_cell
resp(icell) = sum(sum(stas{icell}(:,:,1,itimeSta).*frame));
end

end