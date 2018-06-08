function process_images
% This fcn processes images in the table by applying a flat-field and
% distortion correction on each image. There is no stitching involved in
% this processing procedure.

% Get handles structure
handles=guidata(findall(0,'tag','Nikon_stitch'));

f77=my_fig(77);
axis(f77.s1,'image');
colormap(f77.s1,'bone');
xylabels(f77.s1,'x (\mum)','y (\mum)');
center_axes(f77.s1,'margins',10);

w2=str2double(handles.post_bin.String);
offset=str2double(handles.edit1.String);%intensity offset
T1=str2double(handles.post_thresh.String);%threshold value from IH map
handles.din.IH.ii1=find(handles.din.IH.IH(:)<T1);%pxs indices defined by T1
I_table=handles.uitable1.Data;%extract image files
ii=find(cell2mat(I_table(:,end))==true);%index of images to process
pathnames=I_table(ii,1);%extract out pathnames of imported images
filenames=I_table(ii,2);%extract out filenames of imported images
n=size(filenames,1);%number of files

for dum=1:n
    I=import_tiff_stack([pathnames{dum},filenames{dum}],1,...
        'skip',1,'silence',1);
    
    % Perform offset and flatfield correction
    plane_initial=(I.tiff_stack-offset)./handles.din.IH.IH;
    plane_initial(handles.din.IH.ii1)=0;
    
    % Convert px coord to cartesian coord
    res=I.info.UnknownTags(2).Value;%resolution um/px
    [x,y]=meshgrid(1:size(plane_initial,2),1:size(plane_initial,1));
    x=x.*res;%convert to real position in um
    y=y.*res;
    z1=plane_initial(:);%intensity values
    
    % Perform distortion correction
    if flag==1% correct for distortions
        [x,y]=rm_distort(handles,x,y,res);
    end
    
    % Turn 2d array into column array
    [x,y]=prepareCurveData(x,y);
    
    % Performing binning of current image
    [avg,A,B]=coord2image(x,y,z1,w2,'mean');
    ii1=find(avg(:)>0);
    a=[min(A(ii1)) max(A(ii1))]-min(A(ii1));
    b=[min(B(ii1)) max(B(ii1))]-min(B(ii1));
    
    % Show the processed image
    cla(f77.s1);
    imagesc(f77.s1,'xdata',a,'ydata',b,'cdata',avg);
    drawnow;
    pause(0.5);
end
end

%% Useful fcns
% This fcn handles distortion correction
function [x,y]=rm_distort(handles,x0,y0,res)
    dx=handles.din.D.dx.*res;%x dir distortion
    dy=handles.din.D.dy.*res;%y dir distortion
    x=x0-dx;
    y=y0-dy;
end