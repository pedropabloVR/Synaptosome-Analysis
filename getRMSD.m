function RMSD = getRMSD(X,Y)

% Get mean center of mass
Xc = mean(X);
Yc = mean(Y);

% Calculate root mean squared distance
RMSD = sqrt(mean((X - Xc).^2 + (Y - Yc).^2));

end