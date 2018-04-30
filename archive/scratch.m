handles=guidata(findall(0,'tag','Nikon_stitch'));

a1=get(findall(figure(13),'type','axes'));

img1=handles.prev_plane;
img1b=img1;
img2=handles.cI;

C=str2double(handles.bin_width.String);

[X,Y]=meshgrid(1:size(img1,2),1:size(img1,1));
X1=X.*C+a1.XLim(1)+handles.xpos.I1;
Y1=Y.*C+a1.YLim(1)+handles.ypos.I1;

[X,Y]=meshgrid(1:size(img2,2),1:size(img2,1));
X2=(X).*C+a1.XLim(1)+handles.xpos.I2;
Y2=(Y).*C+a1.YLim(1)+handles.ypos.I2;

nimg=img1-mean(img1(:));
nsec=img2-mean(img1(:));

crr=xcorr2(nimg,nsec);
[ssr,snd] = max(crr(:));
[ij,ji] = ind2sub(size(crr),snd);

[X,Y]=meshgrid(1:size(img2,2),1:size(img2,1));
C=str2double(handles.bin_width.String);
diix=ji-size(img2,2);
diiy=ij-size(img2,1);
X2b=(X+diix);
Y2b=(Y+diiy);


[X,Y]=meshgrid(1:size(img1,2),1:size(img1,1));
X1c=X.*C+a1.XLim(1)+handles.xpos.I1;
Y1c=Y.*C+a1.YLim(1)+handles.ypos.I1;

X2c=X2b.*C+a1.XLim(1)+handles.xpos.I1;
Y2c=Y2b.*C+a1.YLim(1)+handles.ypos.I1;
xpos=(1+diix).*C+a1.XLim(1)+handles.xpos.I1;
ypos=(1+diiy).*C+a1.YLim(1)+handles.ypos.I1;
%% plot images
f1=my_fig(100,{[1 1 1]});
axis(f1.s1,'image');

surf(f1.s1,X1,Y1,img1);
plot(f1.s1,[min(X1(:)) max(X1(:)) max(X1(:)) min(X1(:)) min(X1(:))],...
    [min(Y1(:)) min(Y1(:)) max(Y1(:)) max(Y1(:)) min(Y1(:))],'r-');
surf(f1.s1,X2,Y2,img2);
plot(f1.s1,[min(X2(:)) max(X2(:)) max(X2(:)) min(X2(:)) min(X2(:))],...
    [min(Y2(:)) min(Y2(:)) max(Y2(:)) max(Y2(:)) min(Y2(:))],'g-');

f2=my_fig(101,{[1 1 1]},'marg_h',[0.2 0.1]);
plot(f2.s1,crr(:));
title(f2.s1,'Cross-Correlation');
plot(f2.s1,snd,ssr,'or');

f3=my_fig(102,{[1 1 1]});
axis(f3.s1,'image');
surf(f3.s1,img1);
plot(f3.s1,[1 size(img1,2) size(img1,2) 1 1],...
    [1 1 size(img1,1) size(img1,1) 1],'r-');
surf(f3.s1,X2b,Y2b,img2);

f4=my_fig(103,{[1 1 1]});
axis(f4.s1,'image');
surf(f4.s1,X1c,Y1c,img1);
plot(f4.s1,[min(X1c(:)) max(X1c(:)) max(X1c(:)) min(X1c(:)) min(X1c(:))],...
    [min(Y1c(:)) min(Y1c(:)) max(Y1c(:)) max(Y1c(:)) min(Y1c(:))],'r-');
surf(f4.s1,X2c,Y2c,img2);
plot(f4.s1,[min(X2c(:)) max(X2c(:)) max(X2c(:)) min(X2c(:)) min(X2c(:))],...
    [min(Y2c(:)) min(Y2c(:)) max(Y2c(:)) max(Y2c(:)) min(Y2c(:))],'g-');