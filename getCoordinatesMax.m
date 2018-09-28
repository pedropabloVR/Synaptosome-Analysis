function [idx_x, idx_y] = getCoordinatesMax(X,Y)

idx_y = find(Y == max(Y(:)));
idx_x = X(idx_y);

end