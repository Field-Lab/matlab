function dataa  = compute_nl_for_lnp(WN_datafile,cellIDs,sta_depth,movie_xml,stim_length)
datarun=load_data(WN_datafile);
datarun=load_params(datarun);
datarun=load_neurons(datarun);
datarun=load_sta(datarun);

%cellID = [2146]%cellID_select;
iiicell=1;
for cellID = cellIDs
    cellID
user_STA_depth = sta_depth;
extract_movie_response2;

usta = maskedMovdd*spksGen/sum(spksGen);
input = (usta'*maskedMovdd)';
output = spksGen;


thrr=[];
for iprc=0:2:100
thrr=[thrr;prctile(input,iprc)];
end

fr = [];
in =[];
err=[];
for ipt = 1:length(thrr)-1
    idx = input>thrr(ipt) & input<=thrr(ipt+1);
    fr = [fr;mean(output(idx))];
    err = [err;sqrt(var(output(idx))/sum(idx))];
    in = [in;mean(input(idx))];
end

figure;
errorbar(in,fr,err,'-.');

dataa(iiicell).in=in;
dataa(iiicell).fr=fr;
dataa(iiicell).err=err;
dataa(iiicell).cellID = cellID;
[NN,XX] = hist(input,20);
dataa(iiicell).XX=XX;
dataa(iiicell).NN = NN;
dataa(iiicell).inp_sd = sqrt(var(input));
iiicell=iiicell+1;
end

end