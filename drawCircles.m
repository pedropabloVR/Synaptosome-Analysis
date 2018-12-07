% This function reads in an image and a list of coordinates and draws
% circles with specified radius and colour on the image around the
% coordinates. There are built-in MATLAB functions for this, but if you
% have a really old version of MATLAB or don't have the computer vision
% toolbox, this function is an ok alternative.
% 
% INPUT:
%     img......: rgb image
%     centroids: two-column array with x and y coordinates of centroids
%     radius...: radius of circle
%     colour...: colour triplet (e.g. [0 0 0] for black)
% 
% OUTPUT:
%     img......: rgb image with circles drawn on it
% 
% Based on:
% https://uk.mathworks.com/matlabcentral/fileexchange/39190-draw-a-circle-on-given-image-with-radius
% from Anand Abhishek Singh Modified by Ezra Bruggeman to draw multiple
% circles on an image in any colour.
% 
% Last updated on 16 Sept 2018

function img = drawCircles(img, centroids, radius, colour)

for i = 1:size(centroids,1)
    x = centroids(i,1);
    y = centroids(i,2);
    for j = radius-2:radius+2
        for k = 0:pi/100000:2*pi
            x1 = j*cos(k) + x;
            y1 = j*sin(k) + y;
            p = abs(floor(x1));
            q = abs(floor(y1));
            if((~p==0)&&(~q==0))
                img(p,q,1) = colour(1);
                img(p,q,2) = colour(2);
                img(p,q,3) = colour(3);
            end
        end
    end
end

end