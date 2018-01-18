function yt = smooth_moving_average(y,w,behaviorEdges)
% function yt = smooth_moving_average(y,w,behaviorEdges)
%
% Input variables:
%     "y" is the vector to be smoohted
%     "w" is the half size of the window (default is w=3)
%     "behaviorEdges" the method to use (see below). By defaut the first method
%     will be used.
% The full window has the width 2*w+1 and is centered at the voxel to be
% smoothed.
% There are two possible behaviors at the edges of the vector:
% 1. The vector will be padded with the edge values (default)
% 2. At the edges of the vector, the window is shortened asymetrically

%% allocate variable for output
yt = zeros(size(y));
%% by default use a moving window of 7 sample (i.e. w=3)
ny = length(y);
if ~exist('w','var') || isempty(w)
    w = 3;
end

%% behavior at edges
if ~exist('behaviorEdges', 'var') || isempty(behaviorEdges)
   behaviorEdges = 1; 
end


if behaviorEdges==1
    %% 1. The vector will be padded with the edge values
    if ~isvector(y)
        error('You provided a %dD array as input! Only 1D vectors are allowed!', ndims(y))
    elseif iscolumn(y)
        y = [y(1:w); y; y(end-w+1:end)];
    elseif isrow(y)
        y = [y(1:w) y y(end-w+1:end)];
    end
    for i=w+1:w+ny
        yt(i-w) = nanmean(y(i-w:i+w));
    end
    
else
    %% 2. At the edges of the vector, the window is shortened asymetrically.
    for i=1:length(y)
        k1=i-1;
        k2=length(y)-i;
        if k1>w
            k1 = w;
        end
        if k2>w
            k2 = w;
        end
        yt(i) = nanmean(y(i-k1:i+k2));
    end
    
end