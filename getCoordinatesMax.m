% This function returns the x and y values at the location of the maximum
% of a curve.
% 
% Author: Ezra Bruggeman, Laser Analytics Group
% Last updated on 28 Sept 2018

function [idx_x, idx_y] = getCoordinatesMax(X,Y)

% Get y-value at location of the maximum of the curve
idx_y = find(Y == max(Y(:)));

% Get x-value at location of the maximum of the curve
idx_x = X(idx_y);

end