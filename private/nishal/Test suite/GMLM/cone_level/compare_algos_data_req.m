%% Compare algorithms

%% nnmf 
trainFrac_list =  [0.04,0.06,0.08,0.1,0.2,0.002,0.40,0.004,0.6,0.006,0.8,0.008,0.01,0.012,0.014,0.016,0.018,0.02,1];
trainFrac_list = sort(trainFrac_list,'ascend');


ss_log=[];
tf=[];
ss_log_sh=[];
for itrainfrac = trainFrac_list
try
a=load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data009_nnmf_diff_amt_data/nnmf_3736_su_6_data_%0.02f.mat',itrainfrac));
  tf = [tf;itrainfrac];
ss_log = [ss_log;su_log(4,5)];

mask =ones(size(su_log));
for icone =1:size(mask,1)
mask(icone,icone)=0;
end
mask(4,5)=0;
mask(5,4)=0;
mask = logical(mask);

  ss_log_sh = [ss_log_sh ; max(su_log(mask))];
catch
end

try
load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data009_nnmf_diff_amt_data/nnmf_3736_su_6_data_%0.04f.mat',itrainfrac));
tf = [tf;itrainfrac];
ss_log = [ss_log;su_log(4,5)];

mask =ones(size(su_log));
for icone =1:size(mask,1)
mask(icone,icone)=0;
end
mask(4,5)=0;
mask(5,4)=0;
mask = logical(mask);

  ss_log_sh = [ss_log_sh ; max(su_log(mask))];
catch
end

end
figure;
semilogx(tf,ss_log,'--*')
hold on;
semilogx(tf,ss_log_sh,'--*')
hold on;
%% gmlm
trainFrac_list =  [0.2:-0.05:0.05,0.04:-0.002:0.01,0.6,0.8,1,0.004];
trainFrac_list = sort(trainFrac_list,'ascend');
tf=[];
ss_log=[];
ss_log_sh=[];

mask =ones(size(su_log));
for icone =1:size(mask,1)
mask(icone,icone)=0;
end
mask(4,5)=0;
mask(5,4)=0;
mask = logical(mask);


for itrainfrac = trainFrac_list
try
load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data009_gmlm_diff_amt_data/quad_gmlm_3736_su_6_data_%0.02f.mat',itrainfrac));
    ss_log = [ss_log;su_log(4,5)];
    tf = [tf;itrainfrac];
    ss_log_sh = [ss_log_sh ; max(su_log(mask))];
catch
end

try
load(sprintf('/Volumes/Lab/Users/bhaishahster/Cone_data_alex/2011-12-13-2/data009_gmlm_diff_amt_data/quad_gmlm_3736_su_6_data_%0.05f.mat',itrainfrac));
    ss_log = [ss_log;su_log(4,5)];
    tf = [tf;itrainfrac]; 
    ss_log_sh = [ss_log_sh ; max(su_log(mask))];
catch
end
end

semilogx(tf,ss_log,'--*');
hold on;
semilogx(tf,ss_log_sh,'--*');
xlabel('Fraction of data');
ylabel('Frequency (out of 50)');
legend('NNMF','NNMF other max','GMLM','GMLM other max')
