% Script for testing ad-hoc McStas load routines

% 1D-load:
figure
[x,y,err,xlab,ylab,titl]=Mc1Dload('lmonsrc.sim');
errorbar(x,y,err)
title(titl)
xlabel(xlab)
ylabel(ylab)

% 2D-load:
figure
[x,y,data,z,xlab,ylab,titl]=Mc2Dload('psd.sim');
imagesc(x,y,data);
set(gca,'ydir','normal');
colorbar
xlabel(xlab)
title(titl)
xlabel(xlab)
ylabel(ylab)

% Scanload
figure
[x,y,err,xlab,ylab,titl]=McScanLoad('.',10)
errorbar(x,y,err)
title(titl)
xlabel(xlab)
ylabel(ylab)
