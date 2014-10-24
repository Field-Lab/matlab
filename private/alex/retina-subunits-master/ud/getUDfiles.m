function [histx histy resSUB resSING coneids weights] = getUDfiles(runid,rgcid)

basePath = '/Users/freemanj11/Dropbox/Projects/RetinaSubunits/';


switch runid

	case '2011-10-14-1-data001-0'

		switch rgcid

			case 136

				load(fullfile(basePath,'ud/2011-10-14-1/ud_2011-10-14-1_data003_OFFmidget_f01id136.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-14-1/data001-0/subunit/136-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-14-1/data001-0/subunit/136-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [11 13 19 12];

			case 947

				load(fullfile(basePath,'ud/2011-10-14-1/ud_2011-10-14-1_data003_OFFmidget_f01id947.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-14-1/data001-0/subunit/947-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-14-1/data001-0/subunit/947-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [4 11 10 5];

			case 6811

				load(fullfile(basePath,'ud/2011-10-14-1/ud_2011-10-14-1_data003_OFFmidget_f01id6811.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-14-1/data001-0/subunit/6811-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-14-1/data001-0/subunit/6811-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [12 4 9 1];

			case 7471

				load(fullfile(basePath,'ud/2011-10-14-1/ud_2011-10-14-1_data003_OFFmidget_f01id7471.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-14-1/data001-0/subunit/7471-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-14-1/data001-0/subunit/7471-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [21 10 20 19];

			end

	case '2011-10-25-5-data001-0'

		switch rgcid

			case 184
				load(fullfile(basePath,'ud/2011-10-25-5/ud_2011-10-25-5_data003_OFFmidget_f01id184.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/184-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/184-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [4 3 5 2];

			case 454
				load(fullfile(basePath,'ud/2011-10-25-5/ud_2011-10-25-5_data003_OFFmidget_f01id454.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/454-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/454-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [8 4 5 7];

			case 698
				load(fullfile(basePath,'ud/2011-10-25-5/ud_2011-10-25-5_data003_OFFmidget_f01id698.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/698-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/698-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [4 6 3 1];

			case 948
				load(fullfile(basePath,'ud/2011-10-25-5/ud_2011-10-25-5_data003_OFFmidget_f01id948.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/948-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/948-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [4 6 3 2];

			case 2554
				load(fullfile(basePath,'ud/2011-10-25-5/ud_2011-10-25-5_data003_OFFmidget_f01id2554.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/2554-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/2554-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [7 4 3 6];

			case 2583
				load(fullfile(basePath,'ud/2011-10-25-5/ud_2011-10-25-5_data003_OFFmidget_f01id2583.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/2583-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/2583-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [1 4 7 3];

			case 3169
				load(fullfile(basePath,'ud/2011-10-25-5/ud_2011-10-25-5_data003_OFFmidget_f01id3169.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/3169-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/3169-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [11 3 8 4];

			case 6408
				load(fullfile(basePath,'ud/2011-10-25-5/ud_2011-10-25-5_data003_OFFmidget_f01id6408.mat'));
				clear res
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/6408-subunit-greedy-0.33-fit.mat'));
				resSUB = res;
				load(fullfile(basePath,'anal/2011-10-25-5/data001-0/subunit/6408-subunit-singletons-0.33-fit.mat'));
				resSING = res;
				coneids = [6 5 10 1];

			end


	case '2011-10-25-9-data006-0'

		switch rgcid

		case 4066
			load(fullfile(basePath,'ud/2011-10-25-9/ud_2011-10-25-9_data011_OFFmidget_f06id4066.mat'));
			clear res
			load(fullfile(basePath,'anal/2011-10-25-9/data006-0/subunit/4066-subunit-greedy-0.33-fit.mat'));
			resSUB = res;
			load(fullfile(basePath,'anal/2011-10-25-9/data006-0/subunit/4066-subunit-singletons-0.33-fit.mat'));
			resSING = res;
			coneids = [6 1 3 5];

		case 4489
			load(fullfile(basePath,'ud/2011-10-25-9/ud_2011-10-25-9_data011_OFFmidget_f06id4489.mat'));
			clear res
			load(fullfile(basePath,'anal/2011-10-25-9/data006-0/subunit/4489-subunit-greedy-0.33-fit.mat'));
			resSUB = res;
			load(fullfile(basePath,'anal/2011-10-25-9/data006-0/subunit/4489-subunit-singletons-0.33-fit.mat'));
			resSING = res;
			coneids = [3 10 5 4];

		case 5436
			load(fullfile(basePath,'ud/2011-10-25-9/ud_2011-10-25-9_data011_OFFmidget_f06id5436.mat'));
			clear res
			load(fullfile(basePath,'anal/2011-10-25-9/data006-0/subunit/5436-subunit-greedy-0.33-fit.mat'));
			resSUB = res;
			load(fullfile(basePath,'anal/2011-10-25-9/data006-0/subunit/5436-subunit-singletons-0.33-fit.mat'));
			resSING = res;
			coneids = [1 6 5 4];

		end




	case '2012-04-13-1-data002'

		switch rgcid


		case 2389
			load(fullfile(basePath,'ud/2012-04-13-1/ud_2012-04-13-1_data005_OFFmidget_f02id2389.mat'));
			clear res
			load(fullfile(basePath,'anal/2012-04-13-1/data002/subunit/2536-subunit-greedy-0.33-fit.mat'));
			resSUB = res;
			load(fullfile(basePath,'anal/2012-04-13-1/data002/subunit/2536-subunit-singletons-0.33-fit.mat'));
			resSING = res;
			coneids = [5 1 6 11];

		case 6286
			load(fullfile(basePath,'ud/2012-04-13-1/ud_2012-04-13-1_data005_OFFmidget_f02id6286.mat'));
			clear res
			load(fullfile(basePath,'anal/2012-04-13-1/data002/subunit/6286-subunit-greedy-0.33-fit.mat'));
			resSUB = res;
			load(fullfile(basePath,'anal/2012-04-13-1/data002/subunit/6286-subunit-singletons-0.33-fit.mat'));
			resSING = res;
			coneids = [8 1 4 9];

		case 3002
			load(fullfile(basePath,'ud/2012-04-13-1/ud_2012-04-13-1_data005_OFFmidget_f02id3001.mat'));
			clear res
			load(fullfile(basePath,'anal/2012-04-13-1/data002/subunit/3002-subunit-greedy-0.33-fit.mat'));
			resSUB = res;
			load(fullfile(basePath,'anal/2012-04-13-1/data002/subunit/3002-subunit-singletons-0.33-fit.mat'));
			resSING = res;
			coneids = [4 9 3 1];

		end

	case '2012-09-06-0-data004'

		switch rgcid

		case 7008
			load(fullfile(basePath,'ud/2012-09-06-0/ud_2012-09-06-0_data007_OFFpicks_f04id7008.mat'));
			clear res
			load(fullfile(basePath,'anal/2012-09-06-0/data004/subunit/7008-subunit-greedy-0.33-fit.mat'));
			resSUB = res;
			load(fullfile(basePath,'anal/2012-09-06-0/data004/subunit/7008-subunit-singletons-0.33-fit.mat'));
			resSING = res;
			coneids = [5 6 7 4];

		end

	case '2012-09-13-2-data005'

		switch rgcid

		case 6991
			load(fullfile(basePath,'ud/2012-09-13-2/ud_2012-09-13-2_data007_OFFsparsepicksmidget_f05id6991.mat'));
			clear res
			load(fullfile(basePath,'anal/2012-09-13-2/data001/subunit/6991-subunit-greedy-0.33-fit.mat'));
			resSUB = res;
			load(fullfile(basePath,'anal/2012-09-13-2/data001/subunit/6991-subunit-singletons-0.33-fit.mat'));
			resSING = res;
			coneids = [6 2 5 3];

		end

end

weights.before = rfstimweight.hairyrf;
if isfield(rfstimweight,'stabilitycheck')
	weights.after = rfstimweight.stabilitycheck.hairyrf;
else
	weights.after = [];
end