function varargout = Nikon_stitch(varargin)
% NIKON_STITCH MATLAB code for Nikon_stitch.fig
%      NIKON_STITCH, by itself, creates a new NIKON_STITCH or raises the
%      existing singleton*.
%
%      H = NIKON_STITCH returns the handle to a new NIKON_STITCH or the
%      handle to the existing singleton*.
%
%      NIKON_STITCH('CALLBACK',hObject,eventData,handles,...) calls the
%      local function named CALLBACK in NIKON_STITCH.M with the given input
%      arguments.
%
%      NIKON_STITCH('Property','Value',...) creates a new NIKON_STITCH or
%      raises the existing singleton*.  Starting from the left, property
%      value pairs are applied to the GUI before Nikon_stitch_OpeningFcn
%      gets called.  An unrecognized property name or invalid value makes
%      property application stop.  All inputs are passed to
%      Nikon_stitch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Nikon_stitch

% Last Modified by GUIDE v2.5 08-Jun-2018 09:21:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Nikon_stitch_OpeningFcn, ...
                   'gui_OutputFcn',  @Nikon_stitch_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Nikon_stitch is made visible.
function Nikon_stitch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Nikon_stitch (see VARARGIN)

% Choose default command line output for Nikon_stitch
addpath('format plot fcns/');
addpath('stitching fcns/');

handles.output = hObject;
handles.din.IH_filename=[];
handles.din.pathname=[];
handles.din.IHpathname=[];
handles.uitable2.ColumnName=handles.uitable2.ColumnName(1:5);
handles.uitable2.Data=[];
set(handles.uitable1,'data',[]);
handles.raw_bin_scatter.f=[];
handles.corr_bin_scatter.f=[];
handles.corr_bin_surf.f=[];
handles.corr_no_bin.f=[];
try
    [dev,flag4]=gpu_check();%check for GPU capability
    handles.din.dev=dev;
    handles.din.flag4=flag4;
catch
    disp('Something went wrong with the gpu_check fcn! Stitching will not use GPU.');
    handles.din.dev=[];
    handles.din.flag3=false;
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Nikon_stitch wait for user response (see UIRESUME)
% uiwait(handles.Nikon_stitch);


% --- Outputs from this function are returned to the command line.
function varargout = Nikon_stitch_OutputFcn(hObject, eventdata, handles)  %#ok<*INUSL>
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles) %#ok<*DEFNU>

data=handles.uitable1.Data;

% Get image files
handles.status.String='Select tif images to stitch...';
[filename,pathname,~]=uigetfile('*.tif','Select tif images to stitch',...
    handles.din.pathname,'MultiSelect','on');
handles.din.pathname=pathname;

% If only a single file is collected
if ~iscell(filename)
    filename={filename};
end

row=metadata(pathname,filename);
data=[data;row];
handles.uitable1.Data=data;
handles.Nikon_metadata=row;

guidata(handles.Nikon_stitch,handles);
handles.status.String='Images imported!';

% --- Executes on button press in IH.
function IH_Callback(hObject, eventdata, handles)
% hObject    handle to IH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load flat field map (variable in mat file must be names "IH")
handles.status.String='Select IH correction';
disp('Select IH correction');
if isempty(handles.din.IHpathname)
    def_path=handles.din.pathname;
else
    def_path=handles.din.IHpathname;
end
[filename,pathname,~]=uigetfile('*.mat',...
    'Select calibration image for flat field correction',...
    def_path,'multiselect','off');
handles.din.IHpathname=pathname;
handles.din.IH_filename=filename;
handles.din.IH=load([pathname,filename]);
handles.IH.TooltipString=[pathname,filename];

cla(handles.axes1);
hold(handles.axes1,'on');
contour(handles.axes1,handles.din.IH.IH,'showtext','on','linecolor','k',...
    'fill','on');
axis(handles.axes1,'image');
view(handles.axes1,2);
T=title(handles.axes1,'Flat field correction');
set(T,'units','normalized','position',[0.5 1.3]);
xylabels(handles.axes1,'x (px)','y (px)','fontsize',10);
set(handles.axes1,'dataaspectratio',[1 1 0.001],'fontsize',10);
if handles.rot.Value==1%check if axes need to be rotated
    view(handles.axes1,[90 90]);
    set(handles.axes1,'yaxislocation','right');
else
    view(handles.axes1,[0 90]);
    set(handles.axes1,'yaxislocation','left');
end
select_contour(handles);
guidata(handles.Nikon_stitch,handles);
handles.status.String='IH correction loaded';
disp('IH correction loaded');



% --- Executes on button press in raw.
function raw_Callback(hObject, eventdata, handles)
% hObject    handle to raw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f1=stitch(0,0);
f1.f.Name='Raw stitch';
handles=guidata(handles.Nikon_stitch);
n=size(handles.uitable1.Data,1);%number of files
disp(['Stitching ',num2str(n),' images']);
handles.status.String=['Stitching ',num2str(n),' images'];
for dum=1:n
    In=['I',num2str(dum)];
    coord=handles.coord.(In);%extract out the coordinates
    scatter3(f1.s1,coord(:,1),coord(:,2),coord(:,3),'filled',...
        'cdata',coord(:,3),'sizedata',5,'marker','s');
    
    disp(['stitched image ',num2str(dum)]);
    handles.status.String=['stitched image ',num2str(dum)];
    drawnow;pause(2);
end

% To save on memory, clear handles.coord
handles=rmfield(handles,'coord');

if handles.rot.Value==1
    view(f1.s1,[90 90]);
    set(f1.s1,'yaxislocation','right');
end
set(f1.s1,'clim',[0 2e3],'dataaspectratio',[1 1 10],'zlim',[0 5e3]);
disp('Raw stitching complete');
handles.status.String='Raw stitching complete';
guidata(handles.Nikon_stitch,handles);


% --- Executes on button press in corr.
function corr_Callback(hObject, eventdata, handles)
% hObject    handle to corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

start=clock;

% Perform stitching procedure
f1=stitch(1,0);
f1.f.Name='Flat field and distortion corrected stitch';

handles=guidata(handles.Nikon_stitch);% refresh handle structure

% Create a scatter plot of the stitched images
scatter3(f1.s1,handles.accum_coord(:,1),handles.accum_coord(:,2),...
    handles.accum_coord(:,3),'filled','cdata',handles.accum_coord(:,3),...
    'sizedata',10,'marker','s');

disp('Raw stitching complete');
handles.status.String='Corr stitching complete';
handles.corr_no_bin=f1;%store the figure handle structure
guidata(handles.Nikon_stitch,handles);%update handle structure
rot_Callback(handles.rot,1,handles);%control plotting behavior details

finish=clock;
disp(['Total time: ',num2str(etime(finish,start)),' s']);


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
handles.uitable1.Data=[];
handles.uitable2.Data=[];


% --- Executes on button press in distort.
function distort_Callback(hObject, eventdata, handles)
if isempty(handles.din.IHpathname)
    def_path=handles.din.pathname;
else
    def_path=handles.din.IHpathname;
end

disp('Select distortion file');
handles.status.String='Select distortion file';
%get distortion file
[filename,pathname,~]=uigetfile('*.mat','Select distortion file',...
    def_path);
handles.din.D=load([pathname,filename]);

% Plot the distortion map
[x0,y0]=meshgrid(1:512,1:512);
[x,y]=meshgrid(linspace(1,512,10),linspace(1,512,10));
dxq=interp2(x0,y0,handles.din.D.dx,x,y);
dyq=interp2(x0,y0,handles.din.D.dy,x,y);

quiver(handles.axes2,dxq,dyq);
axis(handles.axes2,'image');
xylabels(handles.axes2,'x (px)','y (px)','fontsize',10);
T=title(handles.axes2,'Distortion');
set(T,'units','normalized','position',[0.5 1.3]);
set(handles.axes2,'fontsize',10);
if handles.rot.Value==1%check if axes need to be rotated
    view(handles.axes2,[90 90]);
    set(handles.axes2,'yaxislocation','right');
else
    view(handles.axes2,[0 90]);
    set(handles.axes2,'yaxislocation','left');
end

guidata(handles.Nikon_stitch,handles);
handles.status.String='Distortion file loaded';
disp('Distortion file loaded');


% --- Executes on button press in bin_stitch.
function bin_stitch_Callback(hObject, eventdata, handles)
start=clock;

% Get metadata from uitable1
I_table=handles.uitable1.Data;%extract image files
ii=cell2mat(I_table(:,end))==true;
filename=I_table(ii,2);%extract out filenames of imported images
n=size(filename,1);%number of files

% Obtain useful user input parameters
w=str2double(handles.post_bin.String);%get the bin size in um
res=handles.uitable1.Data{1,6};%resolution of images

% Update the user on the status
disp(['Stitching ',num2str(n),' images (corrected)']);
handles.status.String=['Stitching ',num2str(n),' images'];

% Perform stitching
f1=stitch(1,1);
f1.f.Name='Binned flat field and distortion corrected stitch';
handles=guidata(handles.Nikon_stitch);%update handles.structure
%                                                                                   Consider using rsetwrite function, so that display performance of the stitched image is better

% Extract coordinates and apply filtering
coord=handles.bin_coord;%xyz coordinates to stitching results
[img,x,y]=coord2image(coord(:,1),coord(:,2),coord(:,3),w,'mean');
img=medfilt2(img,[5 5]);
z=img(:);
ii=find(z<=0);
x(ii)=nan; y(ii)=nan; z(ii)=nan;

% Apply distortion correction to the IH map
select_contour(handles,str2double(handles.post_thresh.String),'corr',1);

% Determine contour lines based on the user specified contour lines drawn
% on the inhomogeneity illumination map
IH=handles.din.IH.IH;
IH=imresize(IH,handles.din.size);
kk=IH(:)>=str2double(handles.post_thresh.String);%logical array
c1=zeros(size(IH));
c1(kk)=1;
c1b=regionprops(c1,'filledimage');
c1=c1b.FilledImage;
c1=imresize(c1,fliplr(handles.din.size));
c2=bwboundaries(c1);
c3=[];
for dum=1:numel(c2)
    c3=cat(1,c3,c2{dum});
end
c3=c3.*(max(x(:))-min(x(:)))/size(img,2);

% Revert the contour line from IH map prior to distortion correction
select_contour(handles,str2double(handles.post_thresh.String));

% Datapoints from the contour line are too clustered, reduce the # indices
% in cc to help with the visibility
freq=1:5:numel(c3(:,1));

W=who;
out_var(W{:});

% Plot 3d scatter plot of the binned results
scatter3(f1.s1,x(:)-min(coord(:,1)),y(:)-min(coord(:,2)),z,'filled',...
'cdata',z,'sizedata',3,'marker','s');

% Create an image plot of the binned stitched results
f2=my_fig(f1.f.Number+10,{[1 1 1]});
set(f2.f,'color','k',...
    'Name','Binned flat field and distortion corrected stitch (image)');

imagesc(f2.s1,'xdata',[min(coord(:,1)) max(coord(:,1))]-min(coord(:,1)),...
    'ydata',[min(coord(:,2)) max(coord(:,2))]-min(coord(:,2)),'cdata',img);

axis(f2.s1,'image');
set(f2.s1,'clim',[0 4095/str2double(handles.post_thresh.String)],...
    'xcolor','w','ycolor','w','zcolor','w','ydir','reverse');
set(f1.s1,'zlim',f2.s1.CLim,'clim',f2.s1.CLim,...
    'xcolor','w','ycolor','w','zcolor','w',...
    'ydir','reverse');%make axis scaling the same
xylabels(f2.s1,'x (\mum)','y (\mum)');

% Store the figure handles in the GUI handle structure
handles.corr_bin_surf=f2;
handles.corr_bin_scatter=f1;

for dum=1:n
    In=['I',num2str(dum)];

    % Show image outlines of final positioning in the stitched image
    plot3(f2.s1,handles.rect.(In)(:,1)-min(coord(:,1)),...
        handles.rect.(In)(:,2)-min(coord(:,2)),...
        ones(size(handles.rect.(In)(:,1))).*f2.s1.ZLim(2),'r-',...
        'linewidth',2);

    plot3(f1.s1,handles.rect.(In)(:,1)-min(coord(:,1)),...
        handles.rect.(In)(:,2)-min(coord(:,2)),...
        ones(size(handles.rect.(In)(:,1))).*f2.s1.ZLim(2),'r-',...
        'linewidth',2);

    % Show outline defined by the thresh of the IH map
    plot3(f2.s1,...
        c3(freq,1)+handles.rect.(In)(1,1)-min(coord(:,1)),...
        c3(freq,2)+handles.rect.(In)(1,2)-min(coord(:,2)),...
        ones(size(c3(freq,1))).*f2.s1.CLim(2),'r.','markersize',4);
    plot3(f1.s1,...
        c3(freq,1)+handles.rect.(In)(1,1)-min(coord(:,1)),...
        c3(freq,2)+handles.rect.(In)(1,2)-min(coord(:,2)),...
        ones(size(c3(freq,1))).*f2.s1.CLim(2),'r.','markersize',4);

    text(f1.s1,min(handles.rect.(In)(:,1))-min(coord(:,1)),...
        min(handles.rect.(In)(:,2))-min(coord(:,2)),...
        f1.s1.ZLim(2),num2str(dum),'color','r',...
        'horizontalalignment','center','verticalalignment','baseline');
    text(f2.s1,min(handles.rect.(In)(:,1))-min(coord(:,1)),...
        min(handles.rect.(In)(:,2))-min(coord(:,2)),...
        f1.s1.ZLim(2),num2str(dum),'color','r',...
        'horizontalalignment','center','verticalalignment','baseline');
end

outline_Callback(handles.outline, 1, handles);
outline_c_Callback(handles.outline_c, 1, handles)
rot_Callback(handles.rot, 1, handles)
ii=find(z>10);
set(f1.s1,'clim',[250 quantile(z(ii),0.95)]);
set(f2.s1,'clim',[250 quantile(z(ii),0.95)]);
    
% Inform user that stitching is complete
guidata(handles.Nikon_stitch,handles);
disp('Corrected bin stitching complete');
handles.status.String='Corr stitching complete';
finish=clock;
disp(['Total time: ',num2str(etime(finish,start)),' s']);


function bin_width_Callback(hObject, eventdata, handles) %#ok<*INUSD>


% --- Executes during object creation, after setting all properties.
function bin_width_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rot. This fcn rotates the axes.
function rot_Callback(hObject, eventdata, handles)

% Find the relevant figures
f1.s1=findall(handles.corr_bin_scatter.f,'type','axes');
f2.s1=findall(handles.corr_bin_surf.f,'type','axes');
f3.s1=findall(handles.corr_no_bin.f,'type','axes');

% Collect all figure handles to be reformatted
all_f=[f1 f2 f3];

for dum=1:numel(all_f)
    
    % if the current axes is an empty graphics holder, skip this interation
    if isempty(all_f(dum).s1)
        continue
    end
    
    colorbar(all_f(dum).s1,'off');
    all_f(dum).c=colorbar(all_f(dum).s1);
    rescale_ax(all_f(dum).s1,hObject.Value);
    RoughCenterAxes(all_f(dum).s1);
    adjust_colorbar(all_f(dum).s1,all_f(dum).c);

    
end


function thresh_Callback(hObject, eventdata, handles)
select_contour(handles,str2double(handles.post_thresh.String));

% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in raw_bin.
function raw_bin_Callback(hObject, eventdata, handles)
% Get metadata from uitable1
I_table=handles.uitable1.Data;%extract image files
ii=find(cell2mat(I_table(:,end))==true);
% w=str2double(handles.bin_width.String);%bin size in microns
pathname=I_table(ii,1);%#ok<NASGU> %extract out pathnames of imported images
filename=I_table(ii,2);%extract out filenames of imported images
n=size(filename,1);%number of files

% Obtain useful user input parameters
% w=str2double(handles.bin_width.String);%get the bin size in um

% Perform stitching of raw images (uncorrected)
f1=stitch(0,2);
f1.f.Name='Binned raw stitch';
handles=guidata(handles.Nikon_stitch);
disp(['Stitching ',num2str(n),' images']);
handles.status.String=['Stitching ',num2str(n),' images'];

%                                                                                   Consider using rsetwrite function, so that display performance of the stitched image is better

% Plot the results (3D scatter plot) of a raw stitched (binned)
coord=handles.raw_bin_coord;%extract out the coordinates
scatter3(f1.s1,coord(:,1),coord(:,2),coord(:,3),'filled',...
    'cdata',coord(:,3),'sizedata',3,'marker','s');

% Format axes
set(f1.s1,'clim',[0 mean(coord(:,3))+2.5.*std(coord(:,3))]);

% Store the figure handles in the GUI handle structure
handles.raw_bin_scatter=f1;

for dum=1:n
    In=['I',num2str(dum)];

    % Show image outlines of final positioning in the stitched image
    plot3(f1.s1,handles.rect.(In)(:,1),handles.rect.(In)(:,2),...
        ones(size(handles.rect.(In)(:,1))).*f1.s1.ZLim(2),'r-','markersize',4);

    % Show outline defined by the thresh of the IH map
    
    text(f1.s1,min(handles.rect.(In)(:,1)),min(handles.rect.(In)(:,2)),...
        f1.s1.ZLim(2),num2str(dum),'color','r',...
        'horizontalalignment','center','verticalalignment','baseline'    );
end

outline_Callback(handles.outline, 1, handles);
rot_Callback(handles.rot, 1, handles)

guidata(handles.Nikon_stitch,handles);
disp('Binned raw stitching complete');
handles.status.String='Raw stitching complete';


% --- Executes on button press in yslice. Obtain results along y direction.
function yslice_Callback(hObject, eventdata, handles)
f1=my_fig(51,{[1 1 1]},'marg_w',[0.2 0.1],'marg_h',[0.2 0.15]);
f1.f.Name='Y slice';
s=str2double(handles.yslice_txt.String);
coord=handles.bin_coord;
ys=unique(coord(:,2));%get unique y coordinates
ys=sort(ys,'ascend');%sort the y coordinates
s2=interp1(ys,ys,s,'nearestneighbor');%find closest coordinate
ii=coord(:,2)==s2;
s_data=coord(ii,:);%extract out the sliced data
scatter3(f1.s1,s_data(:,1),s_data(:,2),s_data(:,3),'filled',...
    'marker','o');
xylabels(f1.s1,'x (\mum)','y (m)','fontweight','bold');
zlabel(f1.s1,'Intensity (a.u.)');
title(f1.s1,['Sliced data at y = ',num2str(s2),' \mum']);
set(findall(f1.s1,'type','text'),'fontweight','bold');
view(f1.s1,[0 0]);

% --- Executes on button press in xslice. Obtain results along x direction.
function xslice_Callback(hObject, eventdata, handles)
f1=my_fig(50,{[1 1 1]},'marg_w',[0.2 0.1],'marg_h',[0.2 0.15]);
f1.f.Name='X slice';
s=str2double(handles.xslice_txt.String);
coord=handles.bin_coord;
xs=unique(coord(:,1));%get unique x coordinates
xs=sort(xs,'ascend');%sort the x coordinates
s2=interp1(xs,xs,s,'nearestneighbor');%find closest coordinate
ii=coord(:,1)==s2;
s_data=coord(ii,:);%extract out the sliced data
scatter3(f1.s1,s_data(:,1),s_data(:,2),s_data(:,3),'filled',...
    'marker','o');
xylabels(f1.s1,'x (\mum)','y (\mum)','fontweight','bold');
zlabel(f1.s1,'Intensity (a.u.)');
title(f1.s1,['Sliced data at x = ',num2str(s2),' \mum']);
set(findall(f1.s1,'type','text'),'fontweight','bold');
view(f1.s1,[90 0]);



% --- Executes on button press in ypos_sort.
function ypos_sort_Callback(hObject, eventdata, handles)
meta=handles.uitable1.Data;
ii=find(cell2mat(meta(:,end))==true);
ii2=find(cell2mat(meta(:,end))==false);
col=cell2mat(meta(ii,4));
if handles.ypos_sort.UserData==1||isempty(handles.ypos_sort.UserData)==1
    [~,ii3]=sort(col,'ascend');
    handles.ypos_sort.UserData=0;
else
    [~,ii3]=sort(col,'descend');
    handles.ypos_sort.UserData=1;
end
handles.uitable1.Data=meta([ii(ii3);ii2],:);
handles.Nikon_metadata=handles.Nikon_metadata([ii(ii3);ii2]);

% --- Executes on button press in xpos_sort.
function xpos_sort_Callback(hObject, eventdata, handles)
meta=handles.uitable1.Data;
ii=find(cell2mat(meta(:,end))==true);
ii2=find(cell2mat(meta(:,end))==false);
col=cell2mat(meta(ii,3));
if handles.xpos_sort.UserData==1||isempty(handles.xpos_sort.UserData)==1
    [~,ii3]=sort(col,'ascend');
    handles.xpos_sort.UserData=0;
else
    [~,ii3]=sort(col,'descend');
    handles.xpos_sort.UserData=1;
end
handles.uitable1.Data=meta([ii(ii3);ii2],:);
handles.Nikon_metadata=handles.Nikon_metadata([ii(ii3);ii2]);


% --- Executes on button press in outline.
function outline_Callback(hObject, eventdata, handles)
try
    if hObject.Value==1
        set(...
            findall(handles.corr_bin_surf.f,'color','r','linestyle','-'),...
            'visible','on');
        set(...
            findall(handles.corr_bin_scatter.f,'color','r','linestyle','-'),...
            'visible','on');
        set(...
            findall(handles.raw_bin_scatter.f,'color','r','linestyle','-'),...
            'visible','on');
        set(...
            findall(handles.corr_bin_surf.f,'color','r','type','text'),...
            'visible','on');
        set(...
            findall(handles.corr_bin_scatter.f,'color','r','type','text'),...
            'visible','on');
        set(...
            findall(handles.raw_bin_scatter.f,'color','r','type','text'),...
            'visible','on');
    else
        set(...
            findall(handles.corr_bin_surf.f,'color','r','linestyle','-'),...
            'visible','off');
        set(...
            findall(handles.corr_bin_scatter.f,'color','r','linestyle','-'),...
            'visible','off');
        set(...
            findall(handles.raw_bin_scatter.f,'color','r','linestyle','-'),...
            'visible','off');
        set(...
            findall(handles.corr_bin_surf.f,'color','r','type','text'),...
            'visible','off');
        set(...
            findall(handles.corr_bin_scatter.f,'color','r','type','text'),...
            'visible','off');
        set(...
            findall(handles.raw_bin_scatter.f,'color','r','type','text'),...
            'visible','off');
    end
catch
    disp('Error with outline_Callback!');
end
rot_Callback(handles.rot, 1, handles)


function TT_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function TT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function post_bin_Callback(hObject, eventdata, handles)
% Notify the user that the plot is refreshing
disp('Refreshing plot with new bin size value');
handles.status.String='Refreshing plot';

try
    f1.f=handles.corr_bin_scatter.f;
    f1.s1=findall(f1.f,'type','axes');
    f2.f=handles.corr_bin_surf.f;
    f2.s1=findall(f2.f,'type','axes');
    cla(f2.s1);
    cla(f1.s1);
catch
    return
end

% Extract out stitched coordinates
X=handles.accum_coord(:,1);
Y=handles.accum_coord(:,2);
Z1=handles.accum_coord(:,3);
w2=str2double(handles.post_bin.String); %extract post_bin value
[avg,A,B]=coord2image(X,Y,Z1,w2,'mean');

%Get rid of zeros
ii1=find(avg(:)>0);
handles.bin_coord=[A(ii1) B(ii1) avg(ii1)];
handles.bin_image=avg;
handles.raw_bin_coord=[A(ii1) B(ii1) avg(ii1)];
handles.raw_bin_image=avg;

%extract out the contour line from IH map
%contour line is drawn after distortion correction
[~,hh]=select_contour(handles,str2double(handles.post_thresh.String),...
    'corr',1);
cc=hh.ContourMatrix';

% Apparently there are spurious datapoints outside the range of the map.
% They need to be removed. Otherwise, the contour outlines will not appear
% correctly in the plots.
iix=find(cc(:,1)<1|cc(:,1)>size(handles.din.IH.IH,2));
iiy=find(cc(:,2)<1|cc(:,2)>size(handles.din.IH.IH,1));
if ~isempty(iix)
    cc(iix,:)=nan;
end
if ~isempty(iiy)
    cc(iiy,:)=nan;
end

% Datapoints from the contour line are too clustered, reduce the # indices
% in cc to help with the visibility
freq=1:10:numel(cc(:,1));

% Revert the contour line from IH map prior to distortion correction
select_contour(handles);

% Plot 3d scatter plot of the binned results
coord=handles.bin_coord;%xyz coordinates to stitching results
scatter3(f1.s1,coord(:,1),coord(:,2),coord(:,3),'filled',...
'cdata',coord(:,3),'sizedata',3,'marker','s');

%                                                                                   For large stitched images, this is likely not a good idea (memory consumption will be large)
% Create a surf plot of the binned stitched results
set(f2.f,'color','k',...
    'Name','Binned flat field and distortion corrected stitch (image)');
[img,~,~]=coord2image(coord(:,1),coord(:,2),coord(:,3),w2,'mean');
imagesc(f2.s1,'xdata',[min(coord(:,1)) max(coord(:,1))],...
    'ydata',[min(coord(:,2)) max(coord(:,2))],'cdata',img);

axis(f2.s1,'image');
set(f2.s1,'clim',[0 4095/str2double(handles.post_thresh.String)],...
    'xcolor','w','ycolor','w','zcolor','w');
set(f1.s1,'zlim',f2.s1.ZLim,'clim',f2.s1.CLim);%make axis scaling the same
xylabels(f2.s1,'x (\mum)','y (\mum)');

% Store the figure handles in the GUI handle structure
handles.corr_bin_surf=f2;
handles.corr_bin_scatter=f1;

for dum=1:size(handles.uitable1.Data,1)
    In=['I',num2str(dum)];
    res=handles.uitable1.Data{dum,6};

    % Show image outlines of final positioning in the stitched image
    plot3(f2.s1,handles.rect.(In)(:,1),handles.rect.(In)(:,2),...
        ones(size(handles.rect.(In)(:,1))).*f2.s1.ZLim(2),'r-');

    plot3(f1.s1,handles.rect.(In)(:,1),handles.rect.(In)(:,2),...
        ones(size(handles.rect.(In)(:,1))).*f2.s1.ZLim(2),'r-');

    % Show outline defined by the thresh of the IH map
    plot3(f2.s1,...
        cc(freq,1).*res+min(handles.rect.(In)(:,1))+w2,...
        cc(freq,2).*res+min(handles.rect.(In)(:,2))+w2,...
        ones(size(cc(freq,1))).*f2.s1.ZLim(2),'r.','markersize',4);
    plot3(f1.s1,...
        cc(freq,1).*res+min(handles.rect.(In)(:,1))+w2,...
        cc(freq,2).*res+min(handles.rect.(In)(:,2))+w2,...
        ones(size(cc(freq,1))).*f2.s1.ZLim(2),'r.','markersize',4);

    text(f1.s1,min(handles.rect.(In)(:,1)),min(handles.rect.(In)(:,2)),...
        f1.s1.ZLim(2),num2str(dum),'color','r',...
        'horizontalalignment','center','verticalalignment','baseline');
    text(f2.s1,min(handles.rect.(In)(:,1)),min(handles.rect.(In)(:,2)),...
        f1.s1.ZLim(2),num2str(dum),'color','r',...
        'horizontalalignment','center','verticalalignment','baseline');
end

outline_Callback(handles.outline, 1, handles);
rot_Callback(handles.rot, 1, handles)

% Inform user that stitching is complete
handles.status.String='Plot refreshed';
disp('Modified (post) bin change complete');
guidata(handles.Nikon_stitch,handles);


% --- Executes during object creation, after setting all properties.
function post_bin_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function post_thresh_Callback(hObject, eventdata, handles)
select_contour(handles,str2double(handles.post_thresh.String));

% --- Executes during object creation, after setting all properties.
function post_thresh_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), ...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nn_sort. Sort rows based on nearest
% neighbor location
function nn_sort_Callback(hObject, eventdata, handles)
raw=handles.uitable1.Data;
x0=cell2mat(raw(:,3));%xpos from uitable1
y0=cell2mat(raw(:,4));%ypos from uitable1
x0=x0-min(x0);
y0=y0-min(y0);
d0=sqrt(x0.^2+y0.^2);

count=1;
index=1;
flag=ones(size(raw,1),1);
dum=1;
while count<numel(flag)
    flag(dum)=nan;
    D=abs(d0-d0(dum)).*flag;
    
    [~,ii]=nanmin(D);
    index=cat(1,index,ii(1));
    dum=ii(1);
    count=count+1;
end
disp(index);
new=raw(index,:);
handles.uitable1.Data=new;


% --- Executes on button press in outline_c.
function outline_c_Callback(hObject, eventdata, handles)
try
    if hObject.Value==1
        set(...
            findall(handles.corr_bin_surf.f,...
            'color','r','linestyle','none'),'visible','on');
        set(...
            findall(handles.corr_bin_scatter.f,...
            'color','r','linestyle','none'),'visible','on');
        set(...
            findall(handles.raw_bin_scatter.f,...
            'color','r','linestyle','none'),'visible','on');
    else
        set(...
            findall(handles.corr_bin_surf.f,...
            'color','r','linestyle','none'),'visible','off');
        set(...
            findall(handles.corr_bin_scatter.f,...
            'color','r','linestyle','none'),'visible','off');
        set(...
            findall(handles.raw_bin_scatter.f,...
            'color','r','linestyle','none'),'visible','off');
    end
catch
    disp('Something is wrong with displaying the outline_c_Callback');
end


% --- Executes on button press in default_pos.
function default_pos_Callback(hObject, eventdata, handles)
I_table=handles.uitable1.Data;%extract image files
pathname=I_table(:,1);%extract out pathnames of imported images
filename=I_table(:,2);%extract out filenames of imported images
tf=I_table(:,end);

% If only a single file is collected
if ~iscell(filename)
    filename={filename};
end
row=metadata(pathname{1},filename);
row(:,end)=tf;
handles.uitable1.Data=row;
handles.Nikon_metadata=row;
disp('Nikon pos restored!');
handles.status.String='Nikon pos restored!';


function out_var(varargin)
% This function output the function variable space to the base workspace
for dum=1:numel(varargin)
    assignin('base',varargin{dum},evalin('caller',varargin{dum}));
end


% --- Executes on button press in step_stitch.
function step_stitch_Callback(hObject, eventdata, handles)

% --- Executes on button press in run_stitch.
function run_stitch_Callback(hObject, eventdata, handles)

% --- Executes on button press in gpu.
function gpu_Callback(hObject, eventdata, handles)
if handles.din.flag4==0
    hObject.Value=0;
    disp('Cannot use GPU accerlation. Prequisites are not met.');
end


% --- Executes on button press in process.
function process_Callback(hObject, eventdata, handles)
process_images;