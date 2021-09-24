
function ambient_geometry_plot(S_R_Geo);

%%%%%%%%%%%%  set ambient source-receiver geometry and plot   %%%%%%%%%%%%%%%%%%%%%
sourceX=S_R_Geo.sx;
sourceY=S_R_Geo.sy;
receiverX=S_R_Geo.rx;
receiverY=S_R_Geo.ry;

[TH_X,R_Y] = cart2pol(sourceX,sourceY);
[TH_S,R_S] = cart2pol(receiverX,receiverY);

figure
polarplot(TH_X, R_Y, 'r*')
hold on;
polarplot(TH_S, R_S, 'bv','MarkerFaceColor','b')


b = findobj(gcf);
c = findall(b,'Type','text');
for phi = 0:30:330
    str = num2str(phi);
    str_new = [str '^\circ'];
    d = findobj(c,'String',str);
    set(d,'String',str_new,'FontSize',12);
end

% figure;
% plot(sourceX, sourceY, 'r*')
% hold on;
% plot(receiverX, receiverY, 'bv','MarkerFaceColor','b')

%
% figure_FontSize=12;
% set(get(gca,'XLabel'),'FontSize',figure_FontSize,'Vertical','top');
% set(get(gca,'YLabel'),'FontSize',figure_FontSize,'Vertical','middle');
% set(findobj('FontSize',10),'FontSize',figure_FontSize);
% set(findobj(get(gca,'Children'),'LineWidth',0.5),'LineWidth',2);
