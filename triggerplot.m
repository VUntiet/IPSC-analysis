function [x,s,t,triggermatrix] = triggerplot(data,triggers,before,after);
 
% function [x,s,t,triggermatrix] = triggerplot(data,triggers,before,after);
% function [x,s,t] = triggerplot(data,triggers,before,after);
% function [x,s] = triggerplot(data,triggers,before,after);
% function [x] = triggerplot(data,triggers,before,after);
% function triggerplot(data,triggers,before,after);
%
% Triggered averaging plot
%
% given 'data' (one dimentional array), the function gives tiggered
% averaging of 'data' at points specified in 'triggers'.  Data points
% considered are 'before' samples before the triggered point and 'after'
% samples after the triggered point (i.e. total of before+after+1=M samples);
% the triggereddata are stored in 'triggermatrix' column (M by triggered
% points matrix), and the average and standard deviation is computed in 'x'
% and 's', respectively.  't' holds a time points, which may be convenient
% for plotting.  
% if no output argument is specified, triggered averaging result with
% the standard error is plotted (see msplot);
%
% 2007.8.10 (debugged minor error - std);



qualified_triggers = sort(triggers(find(triggers>before & triggers<=length(data)-after)));

triggermatrix = zeros(before+after+1,length(qualified_triggers));

for ii=1:length(qualified_triggers)
    triggermatrix(:,ii) = data(qualified_triggers(ii)-before:qualified_triggers(ii)+after);
end

x = mean(triggermatrix,2);
s = std(triggermatrix,0,2);
t = [-before:after];

if nargout==0
    figure;
    msplot(t,triggermatrix);
end
