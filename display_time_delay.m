function varargout = display_time_delay(timerVal)
% function delay = display_time_delay(timerVal)
%
% Function formats the elapsed time for the timer in "timerVal"
% If no timer is provided as input, then the current timer will be used.
% The formated time is either displayed in the command window, or if an output argument is requested, it is returned as output
%
% Input:
%       timerVal - (optional) value of an active timer, as created by fucntion "tic"
%                   If the timerVal argument is not provided, then the current timer will be
%                   used. If none exists, the function will throw an error.
% Output:
%       delay    - (optional) character vector, containing the time delay in the format: 'Elapsed time: hh:mm:ss' (hours:minutes:seconds)
%


%--- retrieve the time delay
if exist('timerVal','var')
    tt = toc(timerVal);
else
    tt = toc;
end
    
%--- split in hours, minutes and seconds
hh=floor(tt/3600);
mm=floor(mod(tt,3600)/60);
ss=round(mod(tt,60));


%--- format and output the time delay
% fprintf('Elapsed time: %dh %dmin %dsec\n', hh, mm, ss)
delay = sprintf('Elapsed time: %02d:%02d:%02d\n', hh, mm, ss);
if nargout>0
    varargout{1} = delay;
else
    fprintf(delay)
end