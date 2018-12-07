% This function estimates the FOV of a SMLM image using the localisations
% and pixelsize. It assumes that there is at least one localisation close
% to the borders of the image.
% 
% Author: Ezra Bruggeman, Laser Analytics Group
% Last updated on 28 Sept 2018

function Fov = estimateFov(X,pixelsize)

X = X/pixelsize;

if X < 64
    Fov = 64;
elseif X < 128
    Fov = 128;
elseif X < 256
    Fov = 256;
else
    Fov = 512;
end

end