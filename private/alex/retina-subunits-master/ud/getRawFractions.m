% 2011-12-13-2
cells1 = {...
'ud_2011-12-13-2_data012-from-data008_data012_data013_data014_OFFmidget_f08id2251',
'ud_2011-12-13-2_data012-from-data008_data012_data013_data014_OFFmidget_f08id3031',
'ud_2011-12-13-2_data012-from-data008_data012_data013_data014_OFFmidget_f08id3586',
'ud_2011-12-13-2_data007-from-data004_data007_data008_OFFMidget_f04id5161',
'ud_2011-12-13-2_data012-from-data008_data012_data013_data014_OFFmidget_f08id7486'};

% 2012-09-24-1
cells2 = {...
'ud_2012-09-24-1_data006_OFFmidget_f03id1861',
'ud_2012-09-24-1_data006_OFFmidget_f03id2807',
'ud_2012-09-24-1_data006_OFFmidget_f03id2837',
'ud_2012-09-24-1_data006_OFFmidget_f03id4186',
'ud_2012-09-24-1_data006_OFFmidget_f03id5761',
'ud_2012-09-24-1_data006_OFFmidget_f03id6826',
'ud_2012-09-24-1_data006_OFFmidget_f03id7397',
'ud_2012-09-24-1_data008_OFFmidget_f05id4126'};


basePath = '/Users/freemanj11/Dropbox/Projects/RetinaSubunits/';

fractionall = [];
for icell = 1:length(cells1)
	load(fullfile(basePath,'ud/2011-12-13-2/',cells1{icell}));
	inds = histx > 0 & histx < 0.3;
	for i=1:4
		for j=1:4
			resp(i,j) = mean(histy{i,j}(inds));
		end
	end
	fraction = [];
	for i=1:4
		for j=i+1:4
			baseline = (resp(i,i) + resp(j,j))/2;
			offresp = (resp(i,j) + resp(j,i))/2;
			fraction = [fraction,offresp/baseline];
		end
	end
	fractionall = [fractionall, fraction];
end

for icell = 1:length(cells2)
	load(fullfile(basePath,'ud/2012-09-24-1/',cells2{icell}));
	inds = histx > 0 & histx < 0.3;
	for i=1:4
		for j=1:4
			resp(i,j) = mean(histy{i,j}(inds));
		end
	end
	fraction = [];
	for i=1:4
		for j=i+1:4
			baseline = (resp(i,i) + resp(j,j))/2;
			offresp = (resp(i,j) + resp(j,i))/2;
			fraction = [fraction,offresp/baseline];
		end
	end
	fractionall = [fractionall, fraction];
end

