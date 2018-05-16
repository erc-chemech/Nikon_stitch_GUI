% This function adjusts the colorbar properly in relation to the axes
function adjust_colorbar(ax,c)
% ax: the axes handle
% c: colorbar handle

xD=diff(ax.XLim);
yD=diff(ax.YLim);

if sum(ax.View==[90 90])==2%if rotated
    r=xD/yD;
else
    r=yD/xD;
end

xD0=ax.Parent.Position(3);
yD0=ax.Parent.Position(4);

colormap(ax,'bone');

% Format colorbar YLabel
c.YLabel.String='Intensity (a.u.)';
c.YLabel.FontSize=ax.YLabel.FontSize;
c.YLabel.Units='pixels';

% Determine positioning of the colorbar
set(ax,'units','pixels');
if r>yD0/xD0
    set(c,'location','eastoutside','color','w',...
        'tickdirection','both','units','pixels',...
        'position',...
        [ax.OuterPosition(1)+ax.OuterPosition(3)+ax.OuterPosition(3)/6,...
        ax.Position(2),20 ax.Position(4)]);
    
    % Determine position of right edge of the colorbar YLabel
    cw=c.YLabel.Extent(1)+c.YLabel.Extent(3);
    
    % Determine the entire plot width
    W=cw+ax.TightInset(3)+ax.TightInset(1)+ax.Position(3)+...
        (c.Position(1)-(ax.Position(1)+ax.Position(3)+ax.TightInset(3)));
    
    
    % Recenter plot area with colorbar
    
    % center based on width
    target=(ax.Parent.InnerPosition(3)-W)/2;
    error_p=ax.Position(1)-ax.TightInset(1)-target;
    ax.Position(1)=ax.Position(1)-error_p;
    c.Position(1)=c.Position(1)-error_p;
    
    % center based on height
    target=(ax.Parent.InnerPosition(4)-...
        (ax.Position(4)+ax.TightInset(2)+ax.TightInset(4)))/2;
    error_p=(ax.Position(2)-ax.TightInset(2))-target;
    ax.Position(2)=ax.Position(2)-error_p;
    c.Position(2)=c.Position(2)-error_p;

else
    set(c,'location','southoutside','color','w',...
        'tickdirection','both','units','pixels',...
        'position',[ax.Position(1),...
        ax.OuterPosition(2)-ax.OuterPosition(2)/6,...
        ax.Position(3) 20]);
    
    % Determine position of bottom edge of the colorbar YLabel
    cw=c.Position(2)+c.YLabel.Extent(2);
    
    % Determine the entire plot height
    W=((c.Position(2)+c.Position(4))-cw)+...
        ax.TightInset(2)+ax.TightInset(4)+ax.Position(4)+...
        ((ax.Position(2)-ax.TightInset(2))-(c.Position(2)+c.Position(4)));
    
    
    % Recenter plot area with colorbar
    
    % center based on height
    target=(ax.Parent.InnerPosition(4)-W)/2;
    error_p=cw-target;
    ax.Position(2)=ax.Position(2)-error_p;
    c.Position(2)=c.Position(2)-error_p;
    
    % center based on width
    target=(ax.Parent.InnerPosition(3)-...
        (ax.Position(3)+ax.TightInset(1)+ax.TightInset(3)))/2;
    error_p=(ax.Position(1)-ax.TightInset(1))-target;
    ax.Position(1)=ax.Position(1)-error_p;
    c.Position(1)=c.Position(1)-error_p;
    
end


set(ax,'units','normalized');
set(c,'units','normalized');
axis(ax,'normal');