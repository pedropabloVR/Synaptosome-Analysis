% This function takes in two images and calculates the Pearson's
% correlation coefficient and Manders' coefficient.
% 
% The PEARSON'S coefficient is sensitive to noise and doesn't give a good
% estimate of colocalization if the stoichiometries are not 1:1 (if the
% signal in the images strongly colocalizes but the average intensities in
% the images are not equal, the Pearson's coefficient will underestimate
% colocalization). The coefficient varies between -1 and 1 (complete
% negative and positive correlation respectively, and 0 standing for no
% correlation).
% The MANDERS' coefficient is based on the Pearson's coefficient, but the
% mean intensities are taken out of the equation. The coefficient varies
% between 0 and 1 (standing for no or complete colocalization).
% The coefficient M1 is defined as the ratio of the 'summed intensities of
% pixels from the first image for which intensity in the other image is
% above zero' to the 'total intensity in the first image'. For M2, the
% images are switched.
%
% INPUT:
%   A, B............: arrays with the same dimensions
% 
% OUTPUT:
%   X.pearsonCoeff..: Pearson's correlation coefficient
%   X.mandersCoeff..: Mander's overlap coefficient
%   X.M1 and X.M2...: Mander's M1 and M2 coefficients
% 
% 
% Author: Ezra Bruggeman, Laser Analytics Group
% Last updated on 22 Sept 2018


function X = getColocCoefficient(A,B)

% Initialize struct
X = struct;

% Pearson's coefficient
ab = (A - mean(A(:))).*(B - mean(B(:)));
aa = (A - mean(A(:))).*(A - mean(A(:)));
bb = (B - mean(B(:))).*(B - mean(B(:)));
X.pearsonCoeff = sum(ab(:))/sqrt(sum(aa(:))*sum(bb(:)));

% Mander's overlap coefficient
ab = A.*B;
aa = A.*A;
bb = B.*B;
X.mandersCoeff = sum(ab(:))/sqrt(sum(aa(:))*sum(bb(:)));
X.M1 = sum(sum(A(B > 0)))/sum(A(:));
X.M2 = sum(sum(B(A > 0)))/sum(B(:));

end