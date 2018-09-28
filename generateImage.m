% Generates a simple SMLM image from a list of localization coordinates.
%
% Draws localizations as points on a black canvas and convolves with
% gaussian kernel with size sigma. The parameter magnification is normally
% equal to 1 and should be ignored. The option to change this value is only
% included so you can visualize very large images faster. But pay
% attention, because coordinates in the resulting image will be scaled by
% this magnification (inverse). For example, visualizing an image that has
% FoV 29960x29960, will become and image with FoV 2996x2996 if
% magnification = 10.
%
% INPUT:
% 	X............: nx1 double, with n the number of localizations
% 	Y............: nx1 double, with n the number of localizations
% 	sigma........: sigma of gaussian kernel localizations are
%                  convolved with
% 	maxInt.......: intensity of localizations before they are
%                  convolved with gaussian filter
% 	Fov..........: size of the image in pixels (e.g. 256, meaning that
%                  image is always a square)
% 	pixelsize....: pixelsize (in nm)
% 	magnification: set to 1, or higher to visualize faster (at the
%                  expense of image quality)
%
% OUTPUT:
%   img..........: nxn uint8 double (with n = Fov*pixelsize/magnification)
% 
% Ideas for improvement:
%   - include a step that pixelates ~ pixelsize
% 
% Author: Ezra Bruggeman, Laser Analytics Group
%
% Last updated on 23 Aug 2018


function img = generateImage(X,Y,sigma,intensity,Fov,pixelsize,magnification)

% Get size of Fov
Fov = ceil((Fov*pixelsize)/magnification) + 1;
%Fov = ceil((Fov*pixelsize)/magnification);

% Create black image with size Fov
img = zeros(Fov,Fov);

% Change pixels in img to a non-zero value (maximum intensity), if they
% contain a localization
for i=1:size(X,1)
    idx_x = ceil(X(i)/magnification);
    idx_y = ceil(Y(i)/magnification);
    if idx_x < 1
        idx_x = 1;
    end
    if idx_y < 1
        idx_y = 1;
    end
    img(idx_x,idx_y) =+ intensity(i);
end

% % Apply gaussian smoothing
% img = imgaussfilt(img,sigma/magnification);

H = fspecial('gaussian',round(3*sigma/magnification),sigma/magnification);
img = imfilter(img,H,'replicate'); 

% Rescale intensities
maxIntensity = max(max(img));
scaleFactor = 256/maxIntensity;
img = uint8(scaleFactor*img);

end