function ecc = getEcc(dataSet)

% ecc = getEcc(dataSet)
%
% lookup table for eccentricities for each piece

switch dataSet
  case '2008-08-27-5/data003'
    ecc = 6.75;
    X = 6.75;
    Y = 0;
  case '2010-03-05-2/data013'
    ecc = 9.5;
    X = 4.75;
    Y = -8.23;
  case '2011-12-13-2/data000-0'
    ecc = 7;
    X = 4.5;
    Y = 4.5;
  case '2011-07-05-2/data002'
    ecc = 5.5;
    X = 5.5;
    Y = 0;
  case '2011-10-25-9/data006-0'
    ecc = 8;
    X = 6.9;
    Y = -4;
  case '2011-10-25-5/data001-0'
    ecc = 8;
    X = 6.9;
    Y = 4;
  case '2011-10-14-1/data001-0'
    ecc = 11;
    X = -5.5;
    Y = -9.5;
  case '2012-04-13-1/data002'
    ecc = 11;
    X = -9.5;
    Y = 5.5;
  case '2011-06-24-6/data005'
    ecc = 8;
    X = 6.9;
    Y = 4;
  case '2011-06-30-0/data003'
    ecc = 8;
    X = 6.9;
    Y = -4;
  case '2012-07-26-1/data004'
    ecc = 5.5;
    X = 5.5;
    Y = 0;
  case '2012-08-09-1/data001'
    ecc = 10.5;
    X = -9.1;
    Y = 5.3;
  case '2012-08-21-0/data001'
    ecc = 11;
    X = -9.5;
    Y = -5.5;
  case '2012-09-06-0/data004'
    ecc = 9; % missing
    X = 9;
    Y = 0;
  case '2012-09-13-2/data001'
    ecc = 5.5;
    X = 5.5;
    Y = 0;
  case '2012-09-18-3/data003'
    ecc = 15;
    X = -15;
    Y = 0;
  case '2012-09-21-2/data007'
    ecc = 15;
    X = -7.5;
    Y = 13;
  case '2012-09-24-1/data003'
    ecc = 11;
    X = 11;
    Y = 0;
end

if X < 0
  ecc = sqrt((X*0.61)^2+Y^2);
else
  ecc = sqrt(X^2+Y^2);
end


%dat = loadData(computer,'2011-06-24-6/data005');
%dat = loadData(computer,'2012-07-26-1/data004');
%dat = loadData(computer,'2011-12-13-2/data000-0');
%dat = loadData(computer,'2012-04-13-1/data002');
%dat = loadData(computer,'2011-10-14-1/data001-0');
%dat = loadData(computer,'2011-10-25-5/data001-0');
%dat = loadData(computer,'2011-10-25-9/data006-0');
%dat = loadData(computer,'2011-12-13-2/data000-0');
%dat = loadData(computer,'vornoi-binary-10-26-2011');
%dat = loadData(computer,'2011-07-05-2/data002');
%dat = loadData(computer,'2011-09-14-1/data000'); % rat
%dat = loadData(computer,'2012-02-09-0/data000'); % rat
%dat = loadData(computer,'2010-09-24-1/data005');
%dat = loadData(computer,'2011-10-25-5/data001-0');
%dat = loadData(computer,'2012-04-13-1/data002');
%dat = loadData(computer,'2011-06-30-0/data003'); % pineapple
%dat = loadData(computer,'2008-08-27-5/data003'); % plantain
%dat = loadData(computer,'2010-03-05-2/data013'); % apple
%dat = loadData(computer,'2008-08-26-2/data001'); % blueberry