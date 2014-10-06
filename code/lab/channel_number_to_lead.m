function leads = channel_number_to_lead(channels);
% CHANNEL_NUMBER_TO_LEAD   Correct channel numbers in 64 electrode data
%
% usage:  leads = channel_number_to_lead(channels);
%
% Correct the channel numbers per LISP code of the same name
% in array.lisp
%
% shlens 2006-06-12
%

% constants from array.lisp
MMM_channel_limit = 64;
RRS_channel_offset = 64;

MMM_channel_leads = [-3 -2 -1 2 13 23 34 45 4 36 52 63 ...
		    10 20 31 42 15 47 53 64 11 21 32 43 ...
		    22 54 51 62 8 19 30 40 41 58 50 61 7 ...
		    18 29 39 26 1 55 60 6 17 28 38 33 44 ...
		    48 59 5 16 27 37 12 49 46 56 3 14 24 35];

for i=1:length(channels)

  if (channels(i) > MMM_channel_limit)
    leads(i) = channels(i) - RRS_channel_offset;
  
  elseif (channels(i) < 0)
    leads(i) = -1 * (abs(channels(i)) - RRS_channel_offset);
  else
    leads(i) = MMM_channel_leads(channels(i)+1);
  end
end