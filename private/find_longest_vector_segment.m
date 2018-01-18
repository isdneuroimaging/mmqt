function [idx1, idx2, lengthSegment] = find_longest_vector_segment(x)
% [idx1, idx2, lengthSegment] = find_longest_vector_segment(x)
%%
x = single(x>0);
t(1,1) = x(1); %--- state
t(1,2) = 1; %--- start index
t(1,3) = nan; %--- stop index
t(1,4) = nan; %--- length
k=1;
for i=2:length(x)
    if x(i)~=t(end,1)
        t(k,3) = i-1;
        t(k,4) = i - t(k,2);
        k=k+1;
        t(k,1) = x(i);
        t(k,2) = i;
    end
end
t(k,3) = length(x);
t(k,4) = length(x) - t(k,2) + 1;

%--- 
tp = t(t(:,1)>0, :);
if ~isempty(tp)
    tpl = tp(tp(:,4)==max(tp(:,4)), : );
    idx1 = tpl(2);
    idx2 = tpl(3);
    lengthSegment = tpl(4);
else
    idx1 = nan;
    idx2 = nan;
    lengthSegment = nan;
end

%---
% plot([idx1 idx1],ylim,'--r');
% plot([idx2 idx2],ylim,'--r');

        
        
        