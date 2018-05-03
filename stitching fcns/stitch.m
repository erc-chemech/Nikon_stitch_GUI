function f1=stitch(flag,flag2)
% Get handles structure
handles=guidata(findall(0,'tag','Nikon_stitch'));

% Flag controls on determining what plots to appear on what figure windows
if flag2==1
    flag3=10;
elseif flag2==2
    flag3=11;
else
    flag3=0;
end

f1=my_fig(flag+flag2+1+flag3);
set(f1.f,'color','k');
set(f1.s1,'xcolor','w','ycolor','w','zcolor','w');
axis(f1.s1,'image');
xylabels(f1.s1,'x (\mum)','y (\mum)');
zlabel(f1.s1,'z (\mum)');

% Do not run this function if the "Run stitch" radio dial is turned off
if handles.run_stitch.Value==0
    return
end

disp([newline,newline,'>>>>>>STITCHING INITIATED<<<<<<']);

handles.accum_coord=[];
Nikon_xpos=handles.Nikon_metadata(:,3);
Nikon_ypos=handles.Nikon_metadata(:,4);
offset=str2double(handles.edit1.String);%intensity offset
T1=str2double(handles.thresh.String);%threshold value from IH map
T2=str2double(handles.post_thresh.String);%threshold value from IH map
handles.din.IH.ii1=find(handles.din.IH.IH(:)<T1);%pxs indices defined by T1
handles.din.IH.ii2=find(handles.din.IH.IH(:)<T2);%pxs indices defined by T2
I_table=handles.uitable1.Data;%extract image files
ii=find(cell2mat(I_table(:,end))==true);
% w=str2double(handles.bin_width.String);%bin size in microns
pathname=I_table(ii,1);%extract out pathnames of imported images
filename=I_table(ii,2);%extract out filenames of imported images
n=size(filename,1);%number of files
for dum=1:n
    In=['I',num2str(dum)];
    I=import_tiff_stack([pathname{dum},filename{dum}],1,...
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
    
    % convert px coord to cartesian coord
    res=I.info.UnknownTags(2).Value;%resolution um/px
    [x,y]=meshgrid(1:size(plane_initial,2),1:size(plane_initial,1));
    x=x.*res;%convert to real position in um
    y=y.*res;
    z1=plane_initial(:);%intensity values
    
    %Perform distortion correction
    if flag==0%do not remove xy distortion
    elseif flag==1%correct for xy distortions
        dx=handles.din.D.dx.*res;%x dir distortion
        dy=handles.din.D.dy.*res;%y dir distortion
        x=x-dx;
        y=y-dy;
    end 
    
    % turn 2d array into column array
    x=x(:);
    y=y(:);
    z1=z1(:);
    
    % Performing binning of current image
    w=ceil(res);
    [cI0,X,Y]=coord2image(x,y,z1,w,'mean');
    
    % Get recorded x and y positions of image
    xposc=handles.uitable1.Data{dum,3};% in um
    yposc=handles.uitable1.Data{dum,4};% in um
    
    % Perform cross correlation for more precise positioning
    if dum>1&&flag==1&&handles.xcorr2.Value==1%do not run on the first iteration
        disp('Performing ZNCC to fine tune positioning...');
        
        % Extract previous image
        pI0=prev.plane;
        
        % Ensure that the current and previous images have the same size
%         size1=[size(pI0,1) size(cI0,1)];
%         size2=[size(pI0,2) size(cI0,2)];
%         pI=zeros([max(size1) max(size2)]);
%         pI(1:size(pI0,1),1:size(pI0,2))=pI0;
%         cI=zeros([max(size1) max(size2)]);
%         cI(1:size(cI0,1),1:size(cI0,2))=cI0;
%         
        % Get recorded x and y positions of the previous image
        xposp=prev.xpos0;
        yposp=prev.ypos0;

        % Determine the relative distance in bin units
        xo=round((xposc-xposp)/w);
        yo=round((yposc-yposp)/w);
%         
%         % Determine overlap regions between the previous and current images
%         if xo>=0&&yo>=0
%             x_overlap=1:(size(cI,2)-xo)-1;
%             y_overlap=1:(size(cI,1)-yo)-1;
%             x_overlap0=(1+xo):size(pI,2);
%             y_overlap0=(1+yo):size(pI,1);
%             ff=1
%         elseif xo>=0&&yo<=0
%             x_overlap=1:(size(cI,2)-xo)-1;
%             y_overlap=(1-yo):size(cI,1);
%             x_overlap0=(1+xo):size(pI,2);
%             y_overlap0=1:(size(pI,1)+yo)-1;
%             ff=2
%         elseif xo<=0&&yo>=0
%             x_overlap=(1-xo):size(cI,2);
%             y_overlap=1:(size(cI,1)-yo)-1;
%             x_overlap0=1:(size(pI,2)+xo)-1;
%             y_overlap0=(1+yo):size(pI,1);
%             ff=3
%         elseif xo<=0&&yo<=0
%             x_overlap=(1-xo):size(cI,2);
%             y_overlap=(1-yo):size(cI,1);
%             x_overlap0=1:(size(pI,2)+xo)-1;
%             y_overlap0=1:(size(pI,1)+yo)-1;
%             ff=4
%         end
%         
%         % Current image overlap area
%         overlap=cI(y_overlap,x_overlap);
%         % Previous image overlap area
%         overlap0=pI(y_overlap0,x_overlap0);
%         
%         W=who;
%         out_var(W{:});
%     
%         % Replace any potential infinity or nan elements with 0
%         ii3a=find(isinf(overlap)|isnan(overlap));
%         ii3b=find(isinf(overlap0)|isnan(overlap0));
%         if ~isempty(ii3a)
%             overlap(ii3a)=0;
%         end
%         if ~isempty(ii3b)
%             overlap0(ii3b)=0;
%         end
%         
%         % Center intensities based on mean value
%         dyn=overlap;
%         static=overlap0;
        



        % Everything is relative to the 1st (previous) image
        outline1=[1 1 size(pI0,2)-1 size(pI0,1)-1];
        
        % Position the current image (2nd image) relative to the
        % coordinates of the previous image (1st image).
        outline2=[1+xo 1+yo size(cI0,2)-1 size(cI0,1)-1];
        
        % Define overlap outline
        a1=max(outline1(1),outline2(1));
        a2=min(outline1(1)+outline1(3),outline2(1)+outline2(3));
        a3=max(outline1(2),outline2(2));
        a4=min(outline1(2)+outline1(4),outline2(2)+outline2(4));
        overlap_outline=[a1 a3 a2-a1 a4-a3];
        
        % Subarea of the previous image (1st image) corresponding to the
        % overlap region
        overlap1=pI0(...
            overlap_outline(2):overlap_outline(2)+overlap_outline(4),...
            overlap_outline(1):overlap_outline(1)+overlap_outline(3));
        
        % Subarea of the current image (2nd iamge) corresponding to the
        % overlap region
        overlap2=cI0(...
            overlap_outline(2)-yo:overlap_outline(2)-yo+overlap_outline(4),...
            overlap_outline(1)-xo:overlap_outline(1)-xo+overlap_outline(3));
        
        % Perform 2D cross correlation (ZNCC)
        static=overlap1;
        dyn=overlap2;
        crrn=xcorr2(dyn,static);
        
        % Determine local region in crrn that is relevant (as defined by
        % TT)
        TT=str2double(handles.TT.String);
        r1=size(dyn,1)-TT;
        r2=size(dyn,1)+TT;
        c1=size(dyn,2)-TT;
        c2=size(dyn,2)+TT;
        %Check to make sure indices are within bounds of crrn size
        if r1<=0; r1=1; end
        if r2>size(crrn,1); r2=size(crrn,1); end
        if c1<=0; c1=1; end
        if c2>size(crrn,2); c2=size(crrn,2); end
        
        crrn2=zeros(size(crrn));
        
        %                                                                           Maybe use column vectors instead and use ind2sub to get peaks
        crrn2(r1:r2,c1:c2)=crrn(r1:r2,c1:c2);
        
        [ypeak, xpeak] = find(crrn2==max(crrn2(:)));%find the maximum
        
        % Calculate the additional offset
        yoffset=-(mean(ypeak)-size(overlap2,1));
        xoffset=-(mean(xpeak)-size(overlap2,2));
        
        % Determine the adjusted xpos and ypos
        xpos=(xo+xoffset)*w+prev.xpos;
        ypos=(yo+yoffset)*w+prev.ypos;
        
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
        
    elseif dum>1&&flag==1&&handles.xcorr2.Value==0
        rel_xpos=xposc-prev.xpos0;
        rel_ypos=yposc-prev.ypos0;
        xpos=prev.xpos+rel_xpos;
        ypos=prev.ypos+rel_ypos;
    elseif dum>1&&flag==0
        xpos=I.info.UnknownTags(5).Value(1);% x pos in um;
        ypos=I.info.UnknownTags(5).Value(2);% y pos in um;
    else
        xpos=xposc;% x pos in um
        ypos=yposc;% y pos in um
    end
    disp(['xpos (\mum): ',num2str(xpos),...
        '   ypos (\mum): ',num2str(ypos)]);
    
    % reimport the original x y z coodinates and apply post threshold value
    % (defined by IH, flat field correction map)
    [x,y]=meshgrid(1:size(plane_ffc,2),1:size(plane_ffc,1));
    x=x.*res;%convert to real position in um
    y=y.*res;
    z1=plane_ffc(:);%intensity values
    
    %Perform distortion correction (again)
    if flag==0%do not remove xy distortion
    elseif flag==1%correct for xy distortions
        dx=handles.din.D.dx.*res;%x dir distortion
        dy=handles.din.D.dy.*res;%y dir distortion
        x=x-dx;
        y=y-dy;
    end 
    
    % Perform stitching
    x=x+xpos;
    y=y+ypos;
    if flag==1%only do this if correction flag is turned on
        z1(handles.din.IH.ii2)=0;
    end
    ii2=find(z1>0);
    x=x(ii2);
    y=y(ii2);
    z1=z1(ii2);
    
    if flag==0||flag2==0
        handles.coord.(In)=[x(:) y(:) z1];
    end
    
    % Update memory of last image stitched
    prev.plane=cI0;
    prev.xpos=xpos;
    prev.ypos=ypos;
    prev.xpos0=xposc;
    prev.ypos0=yposc;
    prev.x=x;
    prev.y=y;
    prev.z1=z1;
    
    
    % Update xpos and ypos on table
    handles.uitable1.Data{dum,3}=xpos;
    handles.uitable1.Data{dum,4}=ypos;
              
    %accumulate all of the scatter3 coordinates
    handles.accum_coord=[handles.accum_coord;x y z1];
    
    % if step stitching is turned one, plot previous and current image
    % stitching
    if handles.step_stitch.Value==1&&dum>1
        In0=['I',num2str(dum-1)];
        cla(handles.axes3);
        scatter3(handles.axes3,handles.accum_coord(:,1),...
            handles.accum_coord(:,2),handles.accum_coord(:,3),...
            'filled','cdata',handles.accum_coord(:,3),...
            'sizedata',3);
        hold(handles.axes3,'on');
        scatter3(handles.axes3,prev.x,...
            prev.y,prev.z1,'filled',...
            'cdata',prev.z1,'sizedata',3);
        view(handles.axes3,2);
        axis(handles.axes3,'image');
    end
    
    disp(['Stitched: ',num2str(dum),' of ',num2str(n)]);
    disp(' ');%Create a blank line in cmd windows  
    
    handles.rect.(In)=...
        [[min(X(:))-w max(X(:))+w max(X(:))+w min(X(:))-w min(X(:))-w]+xpos+w;...
        [min(Y(:))-w min(Y(:))-w max(Y(:))+w max(Y(:))+w min(Y(:))-w]+ypos+w]';
    drawnow;
end

%Performing binning
X=handles.accum_coord(:,1);
Y=handles.accum_coord(:,2);
Z1=handles.accum_coord(:,3);
w2=str2double(handles.post_bin.String);
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
guidata(handles.Nikon_stitch,handles);