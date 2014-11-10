piece = '2012-09-13-2';
run = 'data009';

% define data path
datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarun=load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = set_polarities(datarun);
datarun = load_sta(datarun,'load_sta',[]);

datarunB15 = datarun;
datarunB15 = load_cones(datarunB15,'bayes-msf_15');
datarunB15 = make_mosaic_struct(datarunB15);
datarunB15 = get_sta_fits_from_vision(datarunB15);  
datarunB15 = make_voronoi_masks(datarunB15);

datarunB10=datarun;
datarunB10 = load_cones(datarunB10,'bayes-msf_10');
datarunB10 = make_mosaic_struct(datarunB10);
datarunB10 = get_sta_fits_from_vision(datarunB10);  
datarunB10 = make_voronoi_masks(datarunB10);

datarunB5=datarun;
datarunB5 = load_cones(datarunB5,'bayes-msf_5');
datarunB5 = make_mosaic_struct(datarunB5);
datarunB5 = get_sta_fits_from_vision(datarunB5);  
datarunB5 = make_voronoi_masks(datarunB5);


datarunB15_40 = datarun;
datarunB15_40 = load_cones(datarunB15_40,'bayes-msf_15.00-size40');
datarunB15_40 = make_mosaic_struct(datarunB15_40);
datarunB15_40 = get_sta_fits_from_vision(datarunB15_40);  
datarunB15_40 = make_voronoi_masks(datarunB15_40);



datarunB10_40 = datarun;
datarunB10_40 = load_cones(datarunB10_40,'bayes-msf_10.00-size40');
datarunB10_40 = make_mosaic_struct(datarunB10_40);
datarunB10_40 = get_sta_fits_from_vision(datarunB10_40);  
datarunB10_40 = make_voronoi_masks(datarunB10_40);




a=datarunB15.cones.centers;
b=datarunB15_40.cones.centers;

figure
plot(a(:,1),a(:,2),'+r')
hold on
plot(b(:,1),b(:,2),'x')
legend('bayes 15 size 20','bayes 15 size 40')




a=datarunB10.cones.centers;
b=datarunB10_40.cones.centers;

figure
plot(a(:,1),a(:,2),'+r')
hold on
plot(b(:,1),b(:,2),'x')
legend('bayes 10 size 20','bayes 10 size 40')

