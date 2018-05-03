% Calculate the axes position size and rescale based on the clim and ylim
% This funcion is for the case when the axes is rotated 90 deg.
% counterclock wise.
function rescale_ax2(ax)
view(ax,[90 90]);
set(ax,'units','normalized','yaxislocation','right');
ax.Units='pixels';
xD=diff(ax.XLim);
yD=diff(ax.YLim);
ax.Position(3)=ax.Position(4).*yD./xD;
ax.Position(4)=ax.Position(3).*xD./yD;
ax.Units='normalized';
ax.Position(3:4)=ax.Position(3:4)./max(ax.Position(3:4)).*0.7;
ax.Position(2)=0.5-ax.Position(4)/2;
ax.Position(1)=0.55-ax.Position(3)/2;