function [scaled] = scale_between_0_1(data)
[r,c] = size(data);

min_d = min(data);
max_d = max(data);
range = (max_d - min_d);
scaled = (data - repmat(min_d,r,1)) ./ repmat(range,r,1) ;
end
