% Author: Dr Romain Laine, rfl30@cam.ac.uk, Laser Analytics Group 2015-11-03

function [ r_hist, N_hist, K, N_loc_s, D ] = CalcRipleyK( X1, Y1, X2, Y2, Fov, Area, Analysis_window, r_step )
% For single color Ripley K calculation, simply have CalcRipleyK( X, Y, X, Y, Fov, Area, Analysis_window, r_step )
% Area should already be in nm^2
% The FOV is assummed to be a 1x2 vector (with x and y) 

% How much of the analysis window is removed, it must be >=1
AnWin_rel = 1.1;
r_hist = (0:r_step:Analysis_window)';

% Taking the sample within the region for analysis window
x_s1 = X1((X1 > AnWin_rel*Analysis_window) & (X1 < Fov(1)-AnWin_rel*Analysis_window) & (Y1 > AnWin_rel*Analysis_window) & (Y1 < Fov(2)-AnWin_rel*Analysis_window));
y_s1 = Y1((X1 > AnWin_rel*Analysis_window) & (X1 < Fov(1)-AnWin_rel*Analysis_window) & (Y1 > AnWin_rel*Analysis_window) & (Y1 < Fov(2)-AnWin_rel*Analysis_window));
N_loc_s = length(x_s1);

N_hist = zeros(1,length(r_hist)-1);
for j = 1:N_loc_s
    x_centered = X2 - x_s1(j);
    y_centered = Y2 - y_s1(j);
    r = sqrt(x_centered.^2 + y_centered.^2);
    r(r == 0) = []; % this line is necessary when doing 1 color clustering to get rid of the clustering of the coordinates with themselves
    % r(r > Analysis_window) = []; % This line slows it down, it's faster
    % to use the sub selection of r below
    r_sub = r(r <= Analysis_window);
    
    if ~isempty(r_sub) % if r is empty do nothing
        % N_hist = N_hist + histc(r', r_hist);  % see above
        N_hist = N_hist + histcounts(r_sub', r_hist);
    end
    
end

N_hist = [0; N_hist'];  % add the 0 at r = 0 !
D = length(X2)/Area; % Area should already be in nm^2
K = (1/D)*cumsum(N_hist/N_loc_s); % averaging performed here


% disp(['Number of analysis points: ',num2str(N_loc_s)]);
% disp(['FOV density: ',num2str(10^6*D), ' centroids / um^2']);


end

