f1=my_fig(1);
f2=my_fig(2,{[1 2 1] [1 2 2]});
set(f1.s1,'color','w');
colormap(f1.s1,'bone');

% Everything is relative to the 1st (previous) image
imagesc(f1.s1,[1 size(pI0,2)],[1 size(pI0,1)],pI0,'alphadata',0.5);
outline1=[1 1 size(pI0,2)-1 size(pI0,1)-1];
h1=imrect(f1.s1,outline1);
setColor(h1,'r');

% Position the current image (2nd image) relative to the coordinates of the
% previous image (1st image).
imagesc(f1.s1,[1 size(cI0,2)]+xo,[1 size(cI0,1)]+yo,cI0,'alphadata',0.5);
outline2=[1+xo 1+yo size(cI0,2)-1 size(cI0,1)-1];
h2=imrect(f1.s1,outline2);
setColor(h2,'b');

%define overlap outline
a1=max(outline1(1),outline2(1));
a2=min(outline1(1)+outline1(3),outline2(1)+outline2(3));
a3=max(outline1(2),outline2(2));
a4=min(outline1(2)+outline1(4),outline2(2)+outline2(4));
overlap_outline=[a1 a3 a2-a1 a4-a3];
h3=imrect(f1.s1,overlap_outline);
setColor(h3,'g');

% Subarea of the previous image (1st image) corresponding to the overlap
% region
overlap1=pI0(overlap_outline(2):overlap_outline(2)+overlap_outline(4),...
    overlap_outline(1):overlap_outline(1)+overlap_outline(3));
imagesc(f2.s1,overlap1);

% Subarea of the current image (2nd iamge) corresponding to the overlap
% region
overlap2=cI0(overlap_outline(2)-yo:overlap_outline(2)-yo+overlap_outline(4),...
    overlap_outline(1)-xo:overlap_outline(1)-xo+overlap_outline(3));
imagesc(f2.s2,overlap2);



set(f1.s1,'clim',[0 2e3],'box','off');
axis(f1.s1,'image');
axis([f2.s1 f2.s2],'image');
