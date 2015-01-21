function [day, hour, minute, second] = sec2dhms(sec)
%SEC2DHMS  Convert seconds to days, hours, minutes and seconds.
%
%   [HOUR, MINUTE, SECOND] = SEC2HMS(SEC) converts the number of seconds in
%   SEC into hours, minutes and seconds.

day = fix(sec/(24*3600));
sec = sec - 24*3600*day;
   hour   = fix(sec/3600);      % get number of hours
   sec    = sec - 3600*hour;    % remove the hours
   minute = fix(sec/60);        % get number of minutes
   sec    = sec - 60*minute;    % remove the minutes
   second = sec;
