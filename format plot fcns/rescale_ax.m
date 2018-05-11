% Calculate the axes position size and rescale based on the clim and ylim
% This funcion is for the case when the axes is in the upright position.
function rescale_ax(ax,flag)
%% INPUT VARIABLES
% ax: axes handle
% flag: toggle variable that determines which settings to use (0 or 1)
%
%%

% axes formatting for all cases
set(ax,'units','pixels','tickdir','out','linewidth',1,...
    'yaxislocation','left','clipping','on','box','on');
set(ax.Parent,'units','pixels');
xD=diff(ax.XLim);
yD=diff(ax.YLim);
c=0.7;

if flag==1
    set(ax,'xaxislocation','bottom','view',[90 90]);
    
    ax.Position(3)=c*ax.Parent.Position(4).*yD./xD;
    ax.Position(4)=c*ax.Parent.Position(4);
    
    ax.Units='normalized';
    ax.Position(2)=0.2;
    ax.Position(1)=0.47-ax.Position(3)/2;
    
elseif flag==0
    set(ax,'xaxislocation','top','view',[0 90]);
    
    ax.Position(3)=c*ax.Parent.Position(3);
    ax.Position(4)=c*ax.Parent.Position(3).*yD./xD;
    
    ax.Units='normalized';
    ax.Position(2)=0.55-ax.Position(4)/2;
    ax.Position(1)=0.5-ax.Position(3)/2;
    
end

