function RoughCenterAxes(ax)
% rough centering of axes before colorbar creation(another centering
% procedure will be added after the colorbar has been created for the
% axes).

ax.Units='pixels';
ax.OuterPosition(1)=(ax.Parent.InnerPosition(3)-ax.OuterPosition(3))/2;
ax.OuterPosition(2)=0;
height=ax.Position(4)+ax.Position(2);
ax.OuterPosition(2)=(ax.Parent.InnerPosition(4)-height)/2;
