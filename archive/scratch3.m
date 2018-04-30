close all;
onion   = rgb2gray(imread('onion.png'));
peppers = rgb2gray(imread('peppers.png'));
imshowpair(peppers,onion,'montage')
c = normxcorr2(onion,peppers);
figure, surf(c), shading flat, xlabel('x')

[ypeak, xpeak] = find(c==max(c(:)));
yoffSet = ypeak-size(onion,1);
xoffSet = xpeak-size(onion,2);

figure
imagesc(peppers);
imrect(gca, [xoffSet+1, yoffSet+1, size(onion,2), size(onion,1)]);