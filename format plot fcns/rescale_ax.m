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


if flag==1%rotated by 90 degrees
    set(ax,'xaxislocation','bottom','view',[90 90]);
    
    if yD/xD<=1
        ax.Position(3)=c*ax.Parent.InnerPosition(4)*yD./xD;
        ax.Position(4)=c*ax.Parent.InnerPosition(4);
    else
        ax.Position(4)=c*ax.Parent.InnerPosition(3)*xD/yD;
        ax.Position(3)=c*ax.Parent.InnerPosition(3);
    end
     
elseif flag==0
    set(ax,'xaxislocation','bottom','view',[0 90]);
    
    if yD/xD<=1
        ax.Position(3)=c*ax.Parent.InnerPosition(3);
        ax.Position(4)=c*ax.Parent.InnerPosition(3)*yD/xD;
    else
        ax.Position(3)=c*ax.Parent.InnerPosition(4)*xD/yD;
        ax.Position(4)=c*ax.Parent.InnerPosition(4);
    end
    
end

ax.Units='normalized';
ax.Layer='top';
    


