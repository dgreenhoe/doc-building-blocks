%============================================================================
% Daniel J. Greenhoe
% references:
%   https://octave.org/doc/v4.2.1/Image-Processing.html
%   image: https://www.gentex.com/box/44?iframe=true&width=900&height=600
%          https://www.gentex.com/sites/default/files/imagecache/fire_product/GEC3_24WR.png
%============================================================================
pkg load statistics
image_rgb  = imread("GEC3_24WR.png");
image_gray = rgb2gray(image_rgb);
filter = ones(5,5);
%filter = rand(5,5);
filter = filter / norm(filter);
image_lp   = conv2(image_gray, filter, "same");
%image_lp = image_gray;
image_lp = image_lp / max(max(image_lp))*255;
image_lp = round(image_lp);
size(image_lp)
imwrite(image_rgb,"image_rgb.jpg");
imwrite(image_gray,"image_gray.jpg");
imwrite(image_lp,"image_lp.bmp");
[Dx, Dy]   = gradient (image_lp);