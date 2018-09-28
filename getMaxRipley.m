maxH = max(H);  % Find max value over all elements.
indexOfFirstMax = find(H == maxH, 1, 'first');  % Get first element that is the max.
% Get the x and y values at that index.
maxY = H(indexOfFirstMax);
maxX = r_hist(indexOfFirstMax)

plot(r_hist,H)