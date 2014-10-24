function killList = loadKillList(dataSet,cellType)

switch cellType

	case 'off midget'

	switch dataSet

		case '2008-08-27-5/data003'
			killList = [5539, 4861, 5390, 5616, 5282, 4896, 2181, 287, 170, 467, 6692, 5236, 7471];

		case '2010-03-05-2/data013'
			killList = [6496, 2941];

		case '2011-06-30-0/data003'
			killList = [3679, 860, 3271, 648];

		case '2011-10-14-1/data001-0'
			killList = [];

		case '2011-10-25-5/data001-0'
			killList = [4354, 5091, 842, 1881, 3141, 3263, 5632, 6050];

		case '2011-10-25-9/data006-0'
			killList = [7728, 289, 1771, 1426, 3977, 7550];

		case '2011-12-13-2/data000-0'
			killList = [707, 3948, 3636, 6737, 2633, 812];

		case '2012-04-13-1/data002'
			killList = [4566, 4594, 3124];

		case '2012-07-26-1/data004'
			killList = [3588, 1402, 6876, 3983, 1069, 6726, 616, 291, 7561, 5917, 5333, 5467, 4162];

		case '2012-08-09-1/data001'
			killList = [4789, 3858, 3153, 5661, 5449, 5767, 5149, 1666, 2539, 2013, 3753, 2750, 3062, 3211, 3361, 3783, 5449];

		case '2012-09-06-0/data004'
			killList = [3278, 4036, 4712];

		case '2012-09-13-2/data001'
			killList = [7756, 4640, 5163, 5794, 4865, 4866];

		case '2012-09-18-3/data003'
			killList = [4132, 7401];

		case '2012-09-21-2/data007'
			killList = [571, 7666, 4141];

		case '2012-09-24-1/data003'
			killList = [1772];

		end
			
	case 'on midget'

		switch dataSet

		case '2008-08-27-5/data003'
			killList = [62, 351, 257, 7579, 7279, 6649, 5930, 6048, 5026, 5179];

		case '2010-03-05-2/data013'
			killList = [6649];

		case '2011-06-30-0/data003'
			killList = [6556, 7186];

		case '2011-12-13-2/data000-0'
			killList = [638, 706, 1216, 1291];

		end

	end
end
