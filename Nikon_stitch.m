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

% Last Modified by GUIDE v2.5 02-May-2018 14:37:10

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
handles.output = hObject;
handles.din.IH_filename=[];
handles.din.pathname=[];
handles.din.IHpathname=[];
handles.uitable2.ColumnName=handles.uitable2.ColumnName(1:5);
handles.uitable2.Data=[];
set(handles.uitable1,'data',[]);
handles.corr_bin_scatter.f=[];
handles.corr_bin_surf.f=[];
handles.corr_no_bin.f=[];
addpath('format plot fcns/');
addpath('stitching fcns/');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Nikon_stitch wait for user response (see UIRESUME)
% uiwait(handles.Nikon_stitch);


% --- Outputs from this function are returned to the command line.
function varargout = Nikon_stitch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get data from table
data=handles.uitable1.Data;
offset=size(data,1);

% Get image files
handles.status.String='Select tif images to stitch...';
[filename,pathname,~]=uigetfile('*.tif','Select tif images to stitch',...
    handles.din.pathname,'MultiSelect','on');
handles.din.pathname=pathname;

if iscell(filename)%multiple files selected
    for dum=1:length(filename)
        I=import_tiff_stack([pathname,filename{dum}],1,'skip',1);
        
        x=I.info.UnknownTags(5).Value(1);% x pos in um
        y=I.info.UnknownTags(5).Value(2);% y pos in um
        z=I.info.UnknownTags(5).Value(3);% z pos in um
        res=I.info.UnknownTags(2).Value;%res in um/px
        
        %extract image description
        a=I.info.ImageDescription;
        
        [laser,gain,ex,em,zoom]=metadata(a);
        
        row(dum,:)=...
            {pathname,filename{dum},x,y,z,res,gain,laser,zoom,ex,em,true};
    end
else%single file selected
    In=['I',num2str(1+offset)];
    I=import_tiff_stack([pathname,filename],1,'skip',1);
    
    x=I.info.UnknownTags(5).Value(1);% x pos in um
    y=I.info.UnknownTags(5).Value(2);% y pos in um
    z=I.info.UnknownTags(5).Value(3);% z pos in um
    res=I.info.UnknownTags(2).Value;%res in um/px
 
    %extract image description
    a=I.info.ImageDescription;
 
    [laser,gain,ex,em,zoom]=metadata(a);
 
    row={pathname,filename,x,y,z,res,gain,laser,zoom,ex,em,true};
end

data=[data;row];
handles.uitable1.Data=data;
handles.Nikon_metadata=row;

guidata(handles.Nikon_stitch,handles);
handles.status.String='Images imported!';

function [laser,gain,ex,em,zoom]=metadata(a)
%extract metadata
%laser (%)
ii=strfind(a,'Laser Power')+length('Laser Power}: ');
b=a(ii(1):ii(1)+5);
ii2=strfind(b,newline);
laser=str2double(b(1:ii2-1));%power (%)

%gain (%)
ii=strfind(a,'Gain')+length('Gain}: ');
b=a(ii(1):ii(1)+4);
ii2=strfind(b,' ');
gain=str2double(b(1:ii2-1));%gain

%excitation
ii=strfind(a,'ExcitationWavelength')+length('ExcitationWavelength="');
b=a(ii(1):ii(1)+4);
ii2=strfind(b,'"');
ex=str2double(b(1:ii2-1));%excitation wavelength in nm

%emission
ii=strfind(a,'EmissionWavelength')+length('EmissionWavelength="');
b=a(ii(1):ii(1)+5);
ii2=strfind(b,'"');
em=str2double(b(1:ii2-1));%emission wavelength in nm

%nominal magnification
ii=strfind(a,'NominalMagnification')+length('NominalMagnification="');
b=a(ii(1):ii(1)+4);
ii2=strfind(b,'"');
zoom=str2double(b(1:ii2-1));%zoom

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

% Perform stitching procedure
f1=stitch(1,0);
f1.f.Name='Flat field and distortion corrected stitch';

handles=guidata(handles.Nikon_stitch);% refresh handle structure

% Find images that have been selected for the stitching process
ii=find(cell2mat(handles.uitable1.Data(:,12))==true);

% Create a scatter plot of the stitched images
scatter3(f1.s1,handles.accum_coord(:,1),handles.accum_coord(:,2),...
    handles.accum_coord(:,3),'filled','cdata',handles.accum_coord(:,3),...
    'sizedata',10,'marker','s');

disp('Raw stitching complete');
handles.status.String='Corr stitching complete';
handles.corr_no_bin=f1;%store the figure handle structure
guidata(handles.Nikon_stitch,handles);%update handle structure
rot_Callback(handles.rot,1,handles);%control plotting behavior details


% --- Executes on button press in clear.
function clear_Callback(hObject, eventdata, handles)
% hObject    handle to clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.uitable1.Data=[];
handles.uitable2.Data=[];


% --- Executes on button press in distort.
function distort_Callback(hObject, eventdata, handles)
% hObject    handle to distort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
% hObject    handle to bin_stitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tic

% Get metadata from uitable1
I_table=handles.uitable1.Data;%extract image files
ii=find(cell2mat(I_table(:,end))==true);
filename=I_table(ii,2);%extract out filenames of imported images
n=size(filename,1);%number of files

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
select_contour(handles,str2double(handles.post_thresh.String));

% Obtain useful user input parameters
w=str2double(handles.post_bin.String);%get the bin size in um

% Update the user on the status
disp(['Stitching ',num2str(n),' images (corrected)']);
handles.status.String=['Stitching ',num2str(n),' images'];

% Perform stitching
f1=stitch(1,1);
f1.f.Name='Binned flat field and distortion corrected stitch';
handles=guidata(handles.Nikon_stitch);%update handles.structure

%                                                                                   Consider using rsetwrite function, so that display performance of the stitched image is better
% Plot 3d scatter plot of the binned results
coord=handles.bin_coord;%xyz coordinates to stitching results
scatter3(f1.s1,coord(:,1),coord(:,2),coord(:,3),'filled',...
'cdata',coord(:,3),'sizedata',3,'marker','s');

% Create an image plot of the binned stitched results
f2=my_fig(f1.f.Number+10,{[1 1 1]});
set(f2.f,'color','k',...
    'Name','Binned flat field and distortion corrected stitch (image)');
[img,~,~]=coord2image(coord(:,1),coord(:,2),coord(:,3),w,'mean');
imagesc(f2.s1,'xdata',[min(coord(:,1)) max(coord(:,1))],...
    'ydata',[min(coord(:,2)) max(coord(:,2))],'cdata',img);

axis(f2.s1,'image');
set(f2.s1,'clim',[0 4095/str2double(handles.post_thresh.String)],...
    'xcolor','w','ycolor','w','zcolor','w');
set(f1.s1,'zlim',f2.s1.CLim,'clim',f2.s1.CLim,...
    'xcolor','w','ycolor','w','zcolor','w');%make axis scaling the same
xylabels(f2.s1,['x (\mum)'],['y (\mum)']);

% Store the figure handles in the GUI handle structure
handles.corr_bin_surf=f2;
handles.corr_bin_scatter=f1;

for dum=1:n
    In=['I',num2str(dum)];
    res=handles.uitable1.Data{dum,6};

    % Show image outlines of final positioning in the stitched image
    plot3(f2.s1,handles.rect.(In)(:,1),handles.rect.(In)(:,2),...
        ones(size(handles.rect.(In)(:,1))).*f2.s1.ZLim(2),'r-');

    plot3(f1.s1,handles.rect.(In)(:,1),handles.rect.(In)(:,2),...
        ones(size(handles.rect.(In)(:,1))).*f2.s1.ZLim(2),'r-');

    % Show outline defined by the thresh of the IH map
    plot3(f2.s1,...
        cc(freq,1).*res+min(handles.rect.(In)(:,1))+w,...
        cc(freq,2).*res+min(handles.rect.(In)(:,2))+w,...
        ones(size(cc(freq,1))).*f2.s1.CLim(2),'r.','markersize',4);
    plot3(f1.s1,...
        cc(freq,1).*res+min(handles.rect.(In)(:,1))+w,...
        cc(freq,2).*res+min(handles.rect.(In)(:,2))+w,...
        ones(size(cc(freq,1))).*f2.s1.CLim(2),'r.','markersize',4);

    text(f1.s1,min(handles.rect.(In)(:,1)),min(handles.rect.(In)(:,2)),...
        f1.s1.ZLim(2),num2str(dum),'color','r',...
        'horizontalalignment','center','verticalalignment','baseline');
    text(f2.s1,min(handles.rect.(In)(:,1)),min(handles.rect.(In)(:,2)),...
        f1.s1.ZLim(2),num2str(dum),'color','r',...
        'horizontalalignment','center','verticalalignment','baseline');
end

outline_Callback(handles.outline, 1, handles);
outline_c_Callback(handles.outline_c, 1, handles)
rot_Callback(handles.rot, 1, handles)
    
% Inform user that stitching is complete
guidata(handles.Nikon_stitch,handles);
disp('Corrected bin stitching complete');
handles.status.String='Corr stitching complete';
toc


function bin_width_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function bin_width_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
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
    
    if hObject.Value==1
        rescale_ax2(all_f(dum).s1)
    else
        rescale_ax(all_f(dum).s1)
    end
    
    adjust_colorbar(all_f(dum).s1,all_f(dum).c);
end


function thresh_Callback(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresh as text
%        str2double(get(hObject,'String')) returns contents of thresh as a double
select_contour(handles,str2double(handles.post_thresh.String));

% --- Executes during object creation, after setting all properties.
function thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in raw_bin.
function raw_bin_Callback(hObject, eventdata, handles)
% hObject    handle to raw_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get metadata from uitable1
I_table=handles.uitable1.Data;%extract image files
ii=find(cell2mat(I_table(:,end))==true);
w=str2double(handles.bin_width.String);%bin size in microns
pathname=I_table(ii,1);%extract out pathnames of imported images
filename=I_table(ii,2);%extract out filenames of imported images
n=size(filename,1);%number of files

% Obtain useful user input parameters
w=str2double(handles.bin_width.String);%get the bin size in um

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
% hObject    handle to yslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f1=my_fig(51,{[1 1 1]},'marg_w',[0.2 0.1],'marg_h',[0.2 0.15]);
f1.f.Name='Y slice';
s=str2double(handles.yslice_txt.String);
coord=handles.bin_coord;
ys=unique(coord(:,2));%get unique y coordinates
ys=sort(ys,'ascend');%sort the y coordinates
s2=interp1(ys,ys,s,'nearestneighbor');%find closest coordinate
ii=find(coord(:,2)==s2);
s_data=coord(ii,:);%extract out the sliced data
scatter3(f1.s1,s_data(:,1),s_data(:,2),s_data(:,3),'filled',...
    'marker','o');
xylabels(f1.s1,['x (\mum)'],['y (m)'],'fontweight','bold');
zlabel(f1.s1,'Intensity (a.u.)');
title(f1.s1,['Sliced data at y = ',num2str(s2),' \mum']);
set(findall(f1.s1,'type','text'),'fontweight','bold');
view(f1.s1,[0 0]);

% --- Executes on button press in xslice. Obtain results along x direction.
function xslice_Callback(hObject, eventdata, handles)
% hObject    handle to xslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
f1=my_fig(50,{[1 1 1]},'marg_w',[0.2 0.1],'marg_h',[0.2 0.15]);
f1.f.Name='X slice';
s=str2double(handles.xslice_txt.String);
coord=handles.bin_coord;
xs=unique(coord(:,1));%get unique x coordinates
xs=sort(xs,'ascend');%sort the x coordinates
s2=interp1(xs,xs,s,'nearestneighbor');%find closest coordinate
ii=find(coord(:,1)==s2);
s_data=coord(ii,:);%extract out the sliced data
scatter3(f1.s1,s_data(:,1),s_data(:,2),s_data(:,3),'filled',...
    'marker','o');
xylabels(f1.s1,['x (\mum)'],['y (\mum)'],'fontweight','bold');
zlabel(f1.s1,'Intensity (a.u.)');
title(f1.s1,['Sliced data at x = ',num2str(s2),' \mum']);
set(findall(f1.s1,'type','text'),'fontweight','bold');
view(f1.s1,[90 0]);



% --- Executes on button press in ypos_sort.
function ypos_sort_Callback(hObject, eventdata, handles)
% hObject    handle to ypos_sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
% hObject    handle to xpos_sort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
end
rot_Callback(handles.rot, 1, handles)


function TT_Callback(hObject, eventdata, handles)
% hObject    handle to TT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TT as text
%        str2double(get(hObject,'String')) returns contents of TT as a double


% --- Executes during object creation, after setting all properties.
function TT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function post_bin_Callback(hObject, eventdata, handles)
% hObject    handle to post_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of post_bin as text
%        str2double(get(hObject,'String')) returns contents of post_bin as a double
%Performing binning

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
xylabels(f2.s1,['x (\mum)'],['y (\mum)']);

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
% hObject    handle to post_bin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function post_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to post_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of post_thresh as text
%        str2double(get(hObject,'String')) returns contents of post_thresh as a double
select_contour(handles,str2double(handles.post_thresh.String));

% --- Executes during object creation, after setting all properties.
function post_thresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to post_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
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
    kk=find(flag==1);
    D=abs(d0-d0(dum)).*flag;
    
    [~,ii]=nanmin(D);
    index=[index;ii(1)];
    dum=ii(1);
    count=count+1;
end
disp(index);
new=raw(index,:);
handles.uitable1.Data=new;


% --- Executes on button press in outline_c.
function outline_c_Callback(hObject, eventdata, handles)
% hObject    handle to outline_c (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of outline_c
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
end


% --- Executes on button press in default_pos.
function default_pos_Callback(hObject, eventdata, handles)
% hObject    handle to default_pos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I_table=handles.uitable1.Data;%extract image files
pathname=I_table(:,1);%extract out pathnames of imported images
filename=I_table(:,2);%extract out filenames of imported images
tf=I_table(:,end);

for dum=1:length(filename)
    I=import_tiff_stack([pathname{dum},filename{dum}],1,'skip',1);

    x=I.info.UnknownTags(5).Value(1);% x pos in um
    y=I.info.UnknownTags(5).Value(2);% y pos in um
    z=I.info.UnknownTags(5).Value(3);% z pos in um
    res=I.info.UnknownTags(2).Value;%res in um/px

    %extract image description
    a=I.info.ImageDescription;

    [laser,gain,ex,em,zoom]=metadata(a);

    row(dum,:)=...
        {pathname{dum},filename{dum},x,y,z,res,gain,laser,zoom,ex,em,tf{dum}};    
end

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
