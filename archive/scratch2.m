% crrn=normxcorr2(nimg,nsec);
% f200=my_fig(200,{[1 1 1]});
% surf(f200.s1,crrn);
% xlabel(f200.s1,'x');
% axis(f200.s1,'image');
% set(f200.s1,'dataaspectratio',[1 1 1e-3],'ydir','reverse');
% title(f200.s1,'normxcorr2');
% BW=imregionalmax(crrn);
% 
% [~,snd] = max(crrn(:));%find maximum
%         
% [ypeak,xpeak] = ind2sub(size(crrn),snd);%get row and col subscripts
% yoffset=ypeak-size(cI,1);
% xoffset=xpeak-size(cI,2);

% Properties of current image
[Xc,Yc]=meshgrid(1:size(cI,2),1:size(cI,1));

% Properties of the previous image
[Xp,Yp]=meshgrid(1:size(pI,2),1:size(pI,1));

% Stored (original) rel. positioning in bin units
xposc=I.info.UnknownTags(5).Value(1);
yposc=I.info.UnknownTags(5).Value(2);
xposp=handles.prev_xpos0;
yposp=handles.prev_ypos0;
xo=(xposc-xposp)/w;
yo=(yposc-yposp)/w;

% Determine overlap for current image
xlimc=[min(Xc(:)+xo) max(Xc(:)+xo)];
ylimc=[min(Yc(:)+yo) max(Yc(:)+yo)];

% Determine overlap for previous image
xlimp=[min(Xp(:)) max(Xp(:))];
ylimp=[min(Yp(:)) max(Yp(:))];

% Determin overlap region of current image
iix=find(Xc(:)+xo>xlimp(1)&Xc(:)+xo<xlimp(2));
iiy=find(Yc(:)+yo>ylimp(1)&Yc(:)+yo<ylimp(2));
x_overlap=min(Xc(iix)):max(Xc(iix));
y_overlap=min(Yc(iiy)):max(Yc(iiy));
overlap=zeros(size(pI));
overlap=cI(y_overlap,x_overlap);

overlap0=pI(1:size(overlap,1),1:size(overlap,2));
% overlap=cI(150:172,50:140);
% canvas=zeros(size(pI));
% canvas(1:size(overlap,1),1:size(overlap,2))=overlap;

% TT=2000;
% ii_neg1=find(overlap(:)>TT);
% ii_neg2=find(overlap(:)<=TT);
% overlap_neg=overlap;
% overlap_neg(ii_neg1)=0;
% overlap_neg(ii_neg2)=1;
% 
% ii_neg3=find(pI>TT);
% ii_neg4=find(pI<=TT);
% pI_neg=pI;
% pI_neg(ii_neg3)=0;
% pI_neg(ii_neg4)=1;
% % overlap_neg=reshape(overlap_neg,size(overlap));
% 
% f203=my_fig(203,{[1 1 1]});
% surf(f203.s1,overlap_neg);
% surf(f203.s1,pI_neg);
% xlabel(f203,'x');
% axis(f203.s1,'image');
% view(f203.s1,2);

% [X_overlap,Y_overlap]=meshgrid(x_overlap+xo,y_overlap+yo);

dyn=overlap-mean(overlap(:));
static=overlap0-mean(overlap0(:));

crrn=normxcorr2(dyn,static);
[ypeak, xpeak] = find(crrn==max(crrn(:)));
yoffset=(ypeak-size(overlap,1));
xoffset=(xpeak-size(overlap,2));

[overlapX,overlapY]=meshgrid(1:size(overlap,2),1:size(overlap,1));

f200=my_fig(200,{[1 1 1]});
surf(f200.s1,crrn);
xlabel(f200.s1,'x');
axis(f200.s1,'image');
set(f200.s1,'dataaspectratio',[1 1 1e-3],'ydir','normal');
title(f200.s1,'normxcorr2');

f202=my_fig(202,{[1 1 1]});
surf(f202.s1,Xc+xo+xoffset,Yc+yo+yoffset,cI);
surf(f202.s1,pI);
% surf(f202.s1,overlap0);
% surf(f202.s1,overlap);
% surf(f202.s1,overlapX+xoffset,overlapY+yoffset,overlap);
xlabel(f202.s1,'x');
axis(f202.s1,'image');
set(f202.s1,'clim',[0 5000]);
% imrect(f202.s1, [xoffset+1, yoffset+1, size(overlap,2), size(overlap,1)]);
% set(f202.s1,'xlim',[0 140],'ylim',[0 160]);

% f203=my_fig(203,{[1 1 1]});


% crr=xcorr2(nimg,nsec);
% f201=my_fig(201,{[1 1 1]});
% surf(f201.s1,crr);
% xlabel(f201.s1,'x');
% axis(f201.s1,'image');
% set(f201.s1,'dataaspectratio',[1 1 1e8]);
% title(f201.s1,'xcorr2');


function [out,a1,b1]=coord2image(X,Y,Z1,w,type)
% This fcns converts coordinates of intensities into an image array
n1=min(X):w:max(X);
n2=min(Y):w:max(Y);
a1=linspace(min(X),max(X),numel(n1));
b1=linspace(min(Y),max(Y),numel(n2));

ar=interp1(a1,1:numel(a1),X(:),'nearest');
br=interp1(b1,1:numel(b1),Y(:),'nearest');
if strcmp(type,'mean')
    out=accumarray([br ar],Z1(:),[numel(n2) numel(n1)],@(x) mean(x));
elseif strcmp(type,'none')
    out=accumarray([br ar],Z1(:),[numel(n2) numel(n1)]);
end
end