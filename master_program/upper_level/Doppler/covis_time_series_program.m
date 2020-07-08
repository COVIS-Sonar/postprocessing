% This program is used to build multiple plume property time series using 
% function "covis_time_series".
month{1} = 'jan';
month{2} = 'feb';
month{3} = 'mar';
month{4} = 'apr';
month{5} = 'may';
month{6} = 'jun';
month{7} = 'jul';
month{8} = 'aug';
month{9} = 'sep';
month{10} = 'oct';
month{11} = 'nov';
month{12} = 'dec';
year{1} = '2011';
year{2} = '2012';
year{3} = '2013';
%  matlabpool open local 2
%  parfor i = 9:12
%      covis_time_series('2011',month{i},[1:31],20,32);
%  end
%  matlabpool close
%matlabpool open local 4
for i = 3:12
    covis_time_series('2014',month{i},[1:31],20,32);
end
%matlabpool close
%  matlabpool open local 4
%  parfor i = 1:10
%     covis_time_series('2013',month{i},[1:31],20,32);
%  end
%  matlabpool close