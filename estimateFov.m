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