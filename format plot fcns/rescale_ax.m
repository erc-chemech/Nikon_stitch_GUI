% Calculate the axes position size and rescale based on the clim and ylim
% This funcion is for the case when the axes is in the upright position.
function rescale_ax(ax)
view(ax,[0 90]);
set(ax,'units','normalized','yaxislocation','left');
ax.Units='pixels';
xD=diff(ax.XLim);
yD=diff(ax.YLim);
ax.Position(3)=ax.Position(4).*xD./yD;
ax.Position(4)=ax.Position(3).*yD./xD;
ax.Units='normalized';
ax.Position(3:4)=ax.Position(3:4)./max(ax.Position(3:4)).*0.7;
ax.Position(2)=0.55-ax.Position(4)/2;
ax.Position(1)=0.5-ax.Position(3)/2;