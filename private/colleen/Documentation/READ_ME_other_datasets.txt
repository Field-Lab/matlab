%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Instructions for repeating Malcolm’s analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Datasets:
2007-03-27-1, runs 12-19 (I used runs 14-19)
2007-08-24-4, runs 2-11 (I used runs 4-11)

Scripts:
make_figures_malcolm.m
	Can reproduce all of the plots I made.
motion_script_malcolm.m
	Can reproduce all of the analyses I did, including downsampling. 
	Based on motion_script.m by Marvin Thielke.

Necessary helper scripts:
auto_set_params.m
motion_signal.m
nearest_neighbor_distances.m
pop_motion_signal.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Presentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Malcolm_Presentation.pdf contains a summary of what I did.
All the figures in the presentation can be found in the figures directory.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Notebook searches
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I looked through datasets for other moving bar experiments and object motion experiments. Only monkey experiments are included here. NOTE: Do notebook searches on Bertha! Searching on the locally mounted volumes didn’t work (not all notebooks were searched). I’m not sure why.

I searched using the following command: grep -li ‘<search term>’ `ls -p | grep -v /`
This shows only the filename, ignores case, and only searches files.

I’m listing below the results that came up from the following search terms, and whether or not they’re relevant data sets. Hopefully this will save somebody time at some point.

“moving bar”:
2001-03-16	Not relevant
2001-04-06	Maybe relevant
2001-09-21	Not relevant
2001-09-29	RELEVANT
2001-12-13	RELEVANT
2002-03-15	RELEVANT
2002-04-20	RELEVANT
2002-05-03	RELEVANT
2002-05-10	RELEVANT
2002-07-19	RELEVANT
2002-08-28	RELEVANT
2002-11-14	RELEVANT
2002-12-05	salamander
2003-03-19	guinea pig
2003-03-21	guinea pig
2003-03-24	guinea pig
2003-04-02	RELEVANT
2003-06-10	guinea pig
2003-06-12	guinea pig
2003-06-16	guinea pig
2003-06-17	guinea pig
2003-09-23	guinea pig
2004-01-09	mouse
2004-01-14	rat
2005-01-21	Not relevant
2005-03-25	guinea pig
2005-04-04	guinea pig
2005-04-06	Not relevant
2005-04-14	Not relevant
2005-04-19	Not relevant
2005-04-26	RELEVANT (looks promising, but needs to be cleaned up)
2005-05-02	Guinea pig
2005-05-13	mouse
2005-05-24	Not relevant
2005-05-26	Maybe relevant
2005-06-02	Not relevant
2005-07-07	Not relevant
2005-07-26	Not relevant
2005-08-05	RELEVANT
2005-08-08	Not relevant
2005-09-20	RELEVANT
2005-09-27	RELEVANT
2006-02-10	Not relevant
2006-04-26	Not relevant
2006-05-04	RELEVANT
2006-05-05	Not relevant
2006-05-10	mouse
2006-06-07	mouse
2006-11-08	Not relevant
2006-12-27-A	mouse
2007-03-16	RELEVANT
2007-03-27	RELEVANT (used above)
2007-08-21	RELEVANT
2007-08-24	RELEVANT (used above)
2007-09-18	RELEVANT
2007-12-12	Guinea pig
2007-12-14	Guinea pig
2008-11-09.rtf	Rat
2010-05-19.rtf	Rat
2010-05-20.rtf	Rat
2010-07-22.rtf	Rat
2011-11-18.rtf	Rat
2012-01-20.rtf	Rat
2012-06-18.rtf	Rat
2013-02-07.rtf	Guinea pig
2013-04-01.rtf	Rat

“oms”:
1995-07-20	Not relevant
2003-09-19	Not relevant
2004-03-01	Guinea pig
2005-01-18	Guinea pig
2005-04-04	Guinea pig
2005-05-02	Guinea pig
2005-05-04	Guinea pig
2005-07-25 	Not relevant
2005-08-25-A 	Guinea pig
2005-09-24 	Guinea pig
2005-09-27 	RELEVANT
2005-09-29 	Guinea pig
2006-02-10 	RELEVANT
2006-04-19-A 	Guinea pig
2006-04-24 	RELEVANT
2006-04-26 	RELEVANT
2006-09-07-A	No data
2006-09-22-A 	RELEVANT
2007-02-06 	RELEVANT
2007-03-01-A 	RELEVANT
2008-03-18.rtf	Not relevant
2008-06-10-A.rtf	No data