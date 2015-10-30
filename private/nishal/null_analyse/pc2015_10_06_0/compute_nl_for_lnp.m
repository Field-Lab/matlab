function dataa  = compute_nl_for_lnp(WN_datafile,cellIDs,sta_depth,movie_xml,stim_length)
datarun=load_data(WN_datafile);
datarun=load_params(datarun);
datarun=load_neurons(datarun);
datarun=load_sta(datarun);

%cellID = [2146]%cellID_select;
iiicell=1;
for cellID = cellIDs
    
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
% subplot(1,2,1);
errorbar(in,fr,err,'-.');

dataa(iiicell).in=in;
dataa(iiicell).fr=fr;
dataa(iiicell).err=err;
dataa(iiicell).cellID = cellID;
[NN,XX] = hist(input,20);
dataa(iiicell).XX=XX;
dataa(iiicell).NN = NN;
dataa(iiicell).inp_sd = sqrt(var(input));


% fit the input - firing rate curve and get the slopes at origin and tip
frnz = dataa(iiicell).fr >0; 
f=fit(dataa(iiicell).in(frnz),log(dataa(iiicell).fr(frnz)),'poly2');
g = @(x) exp(f(x));
% subplot(1,2,2);
% plot(g,in,fr);
dataa(iiicell).NL_fit.g=g;
meanin = mean(input);
dataa(iiicell).meaninp = meanin;
dataa(iiicell).NL_fit.g0 = g(meanin) * (2*f.p1*meanin + f.p2);
maxin = max(in);
dataa(iiicell).NL_fit.gend = g(maxin)*(2*f.p1*maxin + f.p2);
dataa(iiicell).NL_fit.nl_idx = log(dataa(iiicell).NL_fit.gend/dataa(iiicell).NL_fit.g0);

iiicell=iiicell+1;
cellID
%pause
end

end