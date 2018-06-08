function f1=stitch(flag,flag2)
%% DESCRIPTION
% This fcn contains the code for the main stitching process.
%
%% INPUT VARIABLES
% flag: toggle variable determining whether or not to perform flat-field
%       corrections and distortion corrections (1=true and 0=false)
% 
% flag2: toggle variable determining whether or not to bin the intensity
% scatter points into a 2d array, based on the meshgrid fcn (1=true and
%       0=false)
% 
%% OUTPUT VARIABLES
% f1: the figure handle in which the stitching results will be shown
% 
%%

% Get handles structure
handles=guidata(findall(0,'tag','Nikon_stitch'));

% Flag3 controls on determining what plots to appear on what figure windows
if flag2==1
    flag3=10;
elseif flag2==2
    flag3=11;
else
    flag3=0;
end

% Figure formatting for output figure
f1=my_fig(flag+flag2+1+flag3);
set(f1.f,'color','k');
axis(f1.s1,'image');
set(f1.s1,'xcolor','w','ycolor','w','zcolor','w');
xylabels(f1.s1,'x (\mum)','y (\mum)');
zlabel(f1.s1,'z (\mum)');

% Do not run this function if the "Run stitch" radio dial is turned off
if handles.run_stitch.Value==0
    return
end

disp([newline,newline,'>>>>>>STITCHING INITIATED<<<<<<']);
flag4=handles.din.flag4;%toggle variable descr. state of useability of gpu
if flag4==true&&handles.gpu.Value==1
    disp('GPU acceleration is turned on');
end
handles.accum_coord=[];
w2=str2double(handles.post_bin.String);
offset=str2double(handles.edit1.String);%intensity offset
T1=str2double(handles.thresh.String);%threshold value from IH map
T2=str2double(handles.post_thresh.String);%threshold value from IH map
handles.din.IH.ii1=find(handles.din.IH.IH(:)<T1);%pxs indices defined by T1
handles.din.IH.ii2=find(handles.din.IH.IH(:)<T2);%pxs indices defined by T2
I_table=handles.uitable1.Data;%extract image files
ii=find(cell2mat(I_table(:,end))==true);%index of images to stitch from table
pathname=I_table(ii,1);%extract out pathnames of imported images
filenames=I_table(ii,2);%extract out filenames of imported images
n=size(filenames,1);%number of files

for dum=1:n
    In=['I',num2str(dum)];
    I=import_tiff_stack([pathname{dum},filenames{dum}],1,...
        'skip',1,'silence',1);
    % Perform offset and flatfield correction
    if flag==0%no flat field correction
        plane_initial=I.tiff_stack-offset;
        plane_ffc=plane_initial;
    elseif flag==1%with flat field correction
        plane_initial=(I.tiff_stack-offset)./handles.din.IH.IH;
        plane_ffc=plane_initial;
        plane_initial(handles.din.IH.ii1)=0;
    end
    
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
    w=ceil(res);% bin based on the resolution of the image
    [cI0,~,~]=coord2image(x,y,z1,w,'none');
    
    % Get recorded x and y positions of image (from uitable)
    xposc=handles.uitable1.Data{dum,3};% in um
    yposc=handles.uitable1.Data{dum,4};% in um
    
    
    % for stitching if user wants to apply cross correlation optimization
    if dum>1&&flag==1&&handles.xcorr2.Value==1
        disp('Performing ZNCC to fine tune positioning...');
        
        % Extract previous image
        pI0=prev.plane;

        % Get recorded x and y positions of the previous image
        xposp=prev.xpos0;
        yposp=prev.ypos0;

        % Determine the relative distance in bin units
        xo=round((xposc-xposp)/w);
        yo=round((yposc-yposp)/w);

        % Everything is relative to the 1st (previous) image
        outline1=[1 1 size(pI0,2)-1 size(pI0,1)-1];% (x y width height)
        
        % Position the current image (2nd image) relative to the
        % coordinates of the previous image (1st image).
        outline2=[1+xo 1+yo size(cI0,2)-1 size(cI0,1)-1];% (x y width height)
        
        % Define overlap outline
        overlap_outline12=overlap_outline(outline1,outline2);
        
        % Subarea of the previous image (1st image) corresponding to the
        % overlap region
        overlap1=pI0(...
            overlap_outline12(2):overlap_outline12(2)+overlap_outline12(4),...
            overlap_outline12(1):overlap_outline12(1)+overlap_outline12(3));
        
        % Subarea of the current image (2nd image) corresponding to the
        % overlap region
        overlap2=cI0(...
            overlap_outline12(2)-yo:...
            overlap_outline12(2)-yo+overlap_outline12(4),...
            overlap_outline12(1)-xo:...
            overlap_outline12(1)-xo+overlap_outline12(3));
        
        % Perform 2D zero normalized cross correlation
        if flag4==false||handles.gpu.Value==0
            tic
            crrn=xcorr2(overlap2-mean(overlap2(:)),overlap1-mean(overlap1(:)));
            toc
        elseif flag4==true&&handles.gpu.Value==1
            tic
            crrn=xcorr2(gpuArray(overlap2-mean(overlap2(:))),...
                gpuArray(overlap1-mean(overlap1(:))));
            wait(handles.din.dev);
            crrn=gather(crrn);%import gpuArray back to 2d array on workspace
            toc
        end
        
        % Determine local region in crrn that is relevant (as defined by
        % TT)
        TT=str2double(handles.TT.String);
        r1=size(overlap2,1)-TT;
        r2=size(overlap2,1)+TT;
        c1=size(overlap2,2)-TT;
        c2=size(overlap2,2)+TT;
        %Check to make sure indices are within bounds of crrn size
        if r1<=0; r1=1; end
        if r2>size(crrn,1); r2=size(crrn,1); end
        if c1<=0; c1=1; end
        if c2>size(crrn,2); c2=size(crrn,2); end
        
        crrn2=zeros(size(crrn));
        crrn2(r1:r2,c1:c2)=crrn(r1:r2,c1:c2);
        
        [ypeak, xpeak] = find(crrn2==max(crrn2(:)),1);%find the maximum
        
        % Calculate the additional offset
        yoffset=-(ypeak-size(overlap2,1));
        xoffset=-(xpeak-size(overlap2,2));
        
        % Determine the adjusted xpos and ypos
        xpos=(xo+xoffset)*w+handles.uitable1.Data{dum-1,3};
        ypos=(yo+yoffset)*w+handles.uitable1.Data{dum-1,4};
        
        disp(['Add. adjust. (x,y in bin units): ',...
            num2str(xoffset),', ',num2str(yoffset)]);
        disp(['repos factor (x,y): ',num2str(abs(xoffset/TT)),...
            ', ',num2str(abs(yoffset/TT))]);
        
        %                                                                            Eventually, I think it will be useful to include what to do if rejected
        if abs(xoffset/TT)>1||abs(yoffset/TT)>1
            disp('Rejected xcorr2! Defaulting to initial rel. pos.');
            rel_xpos=xposc-xposp;
            rel_ypos=yposc-yposp;
            xpos=prev.xpos+rel_xpos;
            ypos=prev.ypos+rel_ypos;
            xoffset=0;
            yoffset=0;
        end
        
        % Show adjustment results in uitable2
        handles.uitable2.Data(dum,:)=...
            [xoffset xoffset*w yoffset yoffset*w,...
            abs(yoffset/TT)+abs(xoffset/TT)];
        handles.status.String=['Stitched: ',num2str(dum),' of ',num2str(n),...
            '    Adjustments (x,y in bin units): ',num2str(xoffset),...
            ', ',num2str(yoffset)];  
        drawnow;

    % for stitching using x and y pos stored in the uitable
    elseif dum>1&&flag==1&&handles.xcorr2.Value==0
        rel_xpos=xposc-prev.xpos0;
        rel_ypos=yposc-prev.ypos0;
        xpos=prev.xpos+rel_xpos;
        ypos=prev.ypos+rel_ypos;
        
    % for stitching using raw Nikon reported positions
    elseif dum>1&&flag==0
        xpos=I.info.UnknownTags(5).Value(1);% x pos in um;
        ypos=I.info.UnknownTags(5).Value(2);% y pos in um;
    
    % Do this for the first image
    else
        xpos=xposc;% x pos in um
        ypos=yposc;% y pos in um
    end
    
    disp(['xpos (\mum): ',num2str(xpos),'   ypos (\mum): ',num2str(ypos)]);
    
    % reimport the original x y z coodinates and apply post threshold value
    % (defined by IH, flat field correction map)
    [x,y]=meshgrid(1:size(plane_ffc,2),1:size(plane_ffc,1));
    x=x.*res;%convert to real position in um
    y=y.*res;
    z1=plane_ffc(:);%intensity values
    
    %Perform distortion correction (again)
    if flag==1% correct for distortions and apply IH filter
        [x,y]=rm_distort(handles,x,y,res);
        z1(handles.din.IH.ii2)=0;
    end
    
    % Update x and y coordinates
    x=x+xpos;
    y=y+ypos;
    
    ii=z1>0;
    x=x(ii);
    y=y(ii);
    z1=z1(ii);
    
    if flag==0||flag2==0
        handles.coord.(In)=[x(:) y(:) z1];
    end    
              
    % accumulate all of the scatter3 coordinates
    handles.accum_coord=[handles.accum_coord;x y z1];
    
    % store outline of the corrected image
    rect_x=[min(x(:)) max(x(:)) max(x(:)) min(x(:)) min(x(:))]';
    rect_y=[min(y(:)) min(y(:)) max(y(:)) max(y(:)) min(y(:))]';
    handles.rect.(In)=[rect_x,rect_y];
    handles.din.size=[numel(min(y(:)):w2:max(y(:))),...
        numel(min(x(:)):w2:max(x(:)))];
    
    % if step stitching is turned on, plot previous and current image
    % stitching
    if handles.step_stitch.Value==1&&dum>1
        cla(handles.axes3);
        scatter3(handles.axes3,x(:),y(:),z1(:),...
            'filled','cdata',z1(:),'sizedata',3);
        hold(handles.axes3,'on');
        scatter3(handles.axes3,prev.x,...
            prev.y,prev.z1,'filled',...
            'cdata',prev.z1,'sizedata',3);
        view(handles.axes3,2);
        axis(handles.axes3,'image');
    end
    
    % Update xpos and ypos on table
    handles.uitable1.Data{dum,3}=xpos;
    handles.uitable1.Data{dum,4}=ypos;
    
    % Update memory of last image stitched
    prev.plane=cI0;
    prev.xpos0=xposc;
    prev.ypos0=yposc;
    prev.xpos=xpos;
    prev.ypos=ypos;
    prev.x=x;
    prev.y=y;
    prev.z1=z1;
    
    disp(['Stitched: ',num2str(dum),' of ',num2str(n)]);
    disp(' ');%Create a blank line in cmd windows  
    
    drawnow;
end

%Performing binning
X=handles.accum_coord(:,1);
Y=handles.accum_coord(:,2);
Z1=handles.accum_coord(:,3);
if flag2>=1    
    [avg,A,B]=coord2image(X,Y,Z1,w2,'mean');

    %Get rid of zeros
    ii1=find(avg(:)>0);
    if flag2==1
        handles.bin_coord=[A(ii1) B(ii1) avg(ii1)];
        handles.bin_image=avg;
    elseif flag2==2
        handles.raw_bin_coord=[A(ii1) B(ii1) avg(ii1)];
        handles.raw_bin_image=avg;
    end
elseif flag2==0
    
end

% Update the handles structure
guidata(handles.Nikon_stitch,handles);

end

%% Useful fcns

% This fcn handles distortion correction
function [x,y]=rm_distort(handles,x0,y0,res)
    dx=handles.din.D.dx.*res;%x dir distortion
    dy=handles.din.D.dy.*res;%y dir distortion
    x=x0-dx;
    y=y0-dy;
end

% Determine the overlap outline region between 2 overlapping rectangles.
% outline1 and outline2 must be a 4-element vector (x y width height)
function overlap_outline=overlap_outline(outline1,outline2)
    a1=max(outline1(1),outline2(1));
    a2=min(outline1(1)+outline1(3),outline2(1)+outline2(3));
    a3=max(outline1(2),outline2(2));
    a4=min(outline1(2)+outline1(4),outline2(2)+outline2(4));
    overlap_outline=[a1 a3 a2-a1 a4-a3];
end


