function plot_section(fieldyz,LAT,DEPTH,latpt,depthpt,obspt,ii,climit,cntrs,txt,colors,cbar)
% function plot_section(fieldyz,LAT,DEPTH,latpt,depthpt,obspt,ii,climit,cntrs,txt,colors,cbar)
%  fieldyz: yz (latitude-depth) section of tracer data
%  LAT: latitude
%  DEPTH: depth    
%  latpt: observational points, latitude
%  depthpt: observational points, depth
%  ii: indices of data to plot
%  climit: colorbar limits
%  cntrs: contour values
%  txt: title
%  colors: color scheme
%  cbar: true for colorbar to be plotted    
    
dcntr = (climit(end)-climit(1))./(size(colors,1))
clrctr = climit(1):dcntr:climit(end)

% determine, based on the color scheme and climit where the lined
% contours should be.
figure(101)
clf(101)
hold on
subplot('position',[.1 .25 .8 .5])
contourf(LAT,-DEPTH,fieldyz,clrctr);
hold on
contour(LAT,-DEPTH,fieldyz,cntrs);
shading flat
[c,h]=contour(LAT,-DEPTH,fieldyz,cntrs,'k');
colormap(colors)
caxis([climit(1) climit(2)])

% Add data points
scatter(latpt(ii),depthpt(ii),40,obspt(ii),'s','filled','MarkerEdgeColor','k');

% manually label the contours
clabel(c,h,'manual')
hold on
set(gca,'XTick',-80:20:80)
xlabel('latitude [$^{\circ}$N]','interpreter','latex','fontsize',11)
set(gca,'YTick',-5000:500:0)
set(gca,'Color',[0.3 0.3 0.3])
set(gca,'YTickLabel',{5,'',4,'',3,'',2,'',1,'',0})
ylabel('depth [km]','fontsize',11)
title(txt,'interpreter','latex','fontsize',11)
if cbar
    h= colorbar;
    %set(h,'YTick',-0.5:0.2:0.5);
end
set(0,'DefaultLineLineWidth',2.0)
set(0,'DefaultTextFontSize',11)
set(0,'DefaultAxesLineWidth',1.6)
set(0,'DefaultAxesFontSize',11)

