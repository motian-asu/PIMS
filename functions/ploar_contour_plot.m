function []=ploar_contour_plot(phi,theta,Z,size_dat)
vals=Z;
az=phi;
incl=theta;

% figure
ph = polarscatter(az * pi/180, incl, [], vals, 'filled'); %convert to radians
ph. LineWidth=0.000001;
ph.SizeData=size_dat; %50
colormap(jet(512));
h = colorbar();
grid off
set(gca, 'FontSize', 36);
set(gca, 'FontName', 'Times New Roman');
set(gca, 'FontWeight', 'bold');
set(gcf,'color','w');

pax = gca;
% pax.ThetaDir='clockwise';
pax.ThetaZeroLocation = 'bottom'%'top';

pxa.ThetaAxis.Visible = 'off'
pax.RTick=[]; % Theta labels
% pxa.YTick=[];
% pxa.YTickLabel = [];
% pxa.YGrid = 'off';
% pxa.ThetaLim=[0 360];
% pxa.Rlim=[0 40];
% ylabel(h, 'value')
end
