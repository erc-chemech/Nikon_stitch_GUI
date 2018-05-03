function [cc,hh]=select_contour(handles,L,varargin)
if nargin==1
    L=str2double(handles.thresh.String);
end

%parse varargin
narginchk(1,inf);
params=inputParser;
params.CaseSensitive=false;
params.addParameter('corr',0,@(x) isnumeric(x));
params.parse(varargin{:});

% Prevent repeated contour lines be plotted in the axes
try
    delete(findall(handles.axes1,'userdata','c'));
end

if params.Results.corr==0%IH contour plotting (no distortion correction)
    [cc,hh]=contour(handles.axes1,handles.din.IH.IH,...
        'levellist',L,...
        'linecolor','r','userdata','c','fill','off','linewidth',2);
elseif params.Results.corr==1%(with distortion correction)
    IH=handles.din.IH.IH;
    dx=handles.din.D.dx;%x dir distortion
    dy=handles.din.D.dy;%y dir distortion
    w=str2double(handles.bin_width.String);% bin size in um
    [X,Y]=meshgrid(1:size(IH,2),1:size(IH,1));
    X=X-dx+w;
    Y=Y-dy+w;
    [cc,hh]=contour(handles.axes1,X,Y,IH,...
        'levellist',L,...
        'linecolor','r','userdata','c','fill','off','linewidth',2);
end
drawnow;